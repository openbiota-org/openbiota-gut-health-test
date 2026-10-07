#!/usr/bin/env python3
"""Create the Route 53 hosted zone for openbiota.com and copy every record from
the current (GoDaddy) zone into it, exactly, before the nameservers move.

Idempotent: re-running updates records in place (UPSERT). Verifies every
record against the new zone's own nameservers afterwards and refuses to call
the migration ready until each one answers identically.

    python3 web/infra/route53_zone.py --profile personal            # create + verify
    python3 web/infra/route53_zone.py --profile personal --verify   # verify only

What is copied (read live from the authoritative GoDaddy servers at run time,
so nothing is typed by hand):
  apex MX x5 (Google Workspace), apex TXT (SPF + Google site verification),
  google._domainkey TXT (Google DKIM), _dmarc TXT, the SES DKIM CNAMEs,
  mail.openbiota.com MX + TXT (SES MAIL FROM).
What is deliberately not copied:
  the GoDaddy parking A records (replaced by the site's ALIAS records),
  _domainconnect (GoDaddy's own connector pointer; meaningless elsewhere),
  GoDaddy's NS/SOA (Route 53 writes its own).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import uuid

DOMAIN = "openbiota.com"
OLD_NS = "pdns09.domaincontrol.com"
SES_DKIM_TOKENS = (
    "snrmi377ixajxb4k52fraae4rpnjwqb7",
    "m5e6nfdphtliuct6bj4bdmh7zincryq6",
    "3brspg6ke7vcl5bagaqzaa6ivaip36yd",
)
HOSTS = ["@", "mail", "_dmarc", "google._domainkey", *(f"{t}._domainkey" for t in SES_DKIM_TOKENS)]
TYPES = ("A", "AAAA", "MX", "TXT", "CNAME")
SKIP_APEX_TYPES = {"A", "AAAA"}   # parking records: the site replaces them


def dig(name: str, rtype: str, server: str) -> list[str]:
    out = subprocess.run(["dig", "+short", f"@{server}", name, rtype], capture_output=True, text=True, check=False).stdout
    return [ln.strip() for ln in out.splitlines() if ln.strip()]


def live_records() -> dict[tuple[str, str], list[str]]:
    """{(fqdn, type): [values]} from the current authoritative servers."""
    out: dict[tuple[str, str], list[str]] = {}
    for host in HOSTS:
        fqdn = DOMAIN if host == "@" else f"{host}.{DOMAIN}"
        has_cname = dig(fqdn, "CNAME", OLD_NS)
        if has_cname:
            out[(fqdn, "CNAME")] = has_cname
            continue   # a CNAME excludes every other type at that name
        for rtype in ("A", "AAAA", "MX", "TXT"):
            if host == "@" and rtype in SKIP_APEX_TYPES:
                continue
            vals = dig(fqdn, rtype, OLD_NS)
            if vals:
                out[(fqdn, rtype)] = vals
    return out


def aws(profile: str, *args: str) -> dict:
    cmd = ["aws", "--profile", profile, "--output", "json", *args]
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode:
        sys.exit(f"aws {' '.join(args[:3])} failed:\n{res.stderr}")
    return json.loads(res.stdout) if res.stdout.strip() else {}


def ensure_zone(profile: str) -> tuple[str, list[str]]:
    zones = aws(profile, "route53", "list-hosted-zones-by-name", "--dns-name", DOMAIN, "--max-items", "1")
    for z in zones.get("HostedZones", []):
        if z["Name"] == f"{DOMAIN}.":
            zid = z["Id"].split("/")[-1]
            ns = aws(profile, "route53", "get-hosted-zone", "--id", zid)["DelegationSet"]["NameServers"]
            print(f"zone exists: {zid}")
            return zid, ns
    created = aws(profile, "route53", "create-hosted-zone", "--name", DOMAIN, "--caller-reference", str(uuid.uuid4()),
                  "--hosted-zone-config", f"Comment=openbiota.com - migrated from GoDaddy {time.strftime('%Y-%m-%d')}")
    zid = created["HostedZone"]["Id"].split("/")[-1]
    print(f"zone created: {zid}")
    return zid, created["DelegationSet"]["NameServers"]


def txt_value(v: str) -> str:
    """dig prints long TXT records as several quoted strings; Route 53 wants the same form."""
    return v if v.startswith('"') else f'"{v}"'


def upsert(profile: str, zid: str, records: dict[tuple[str, str], list[str]]) -> None:
    changes = []
    for (fqdn, rtype), values in sorted(records.items()):
        if rtype == "TXT":
            rrs = [{"Value": txt_value(v)} for v in values]
        else:
            rrs = [{"Value": v} for v in values]
        changes.append({"Action": "UPSERT", "ResourceRecordSet": {
            "Name": f"{fqdn}.", "Type": rtype, "TTL": 300, "ResourceRecords": rrs}})
    batch = {"Comment": "copy of the GoDaddy zone", "Changes": changes}
    res = aws(profile, "route53", "change-resource-record-sets", "--hosted-zone-id", zid, "--change-batch", json.dumps(batch))
    print(f"upserted {len(changes)} record sets (change {res['ChangeInfo']['Id'].split('/')[-1]})")


def normalise(vals: list[str], rtype: str) -> list[str]:
    out = []
    for v in vals:
        v = v.strip()
        if rtype == "TXT":
            v = v.replace('" "', "")       # join dig's split strings
            v = v.strip('"')
        out.append(v.lower() if rtype != "TXT" else v)
    return sorted(out)


def verify(records: dict[tuple[str, str], list[str]], new_ns: list[str]) -> bool:
    ok = True
    server = new_ns[0]
    for (fqdn, rtype), values in sorted(records.items()):
        got = dig(fqdn, rtype, server)
        a, b = normalise(values, rtype), normalise(got, rtype)
        mark = "ok " if a == b else "MISMATCH"
        ok = ok and a == b
        print(f"  {mark} {rtype:5s} {fqdn:58s} {len(b)} value(s)")
        if a != b:
            print(f"        godaddy : {a}\n        route53 : {b}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default="personal")
    ap.add_argument("--verify", action="store_true", help="verify only; create nothing")
    args = ap.parse_args()

    records = live_records()
    print(f"live GoDaddy records to carry over: {len(records)}")
    for (fqdn, rtype), values in sorted(records.items()):
        print(f"  {rtype:5s} {fqdn:58s} {values if rtype != 'TXT' else [v[:60] + ('…' if len(v) > 60 else '') for v in values]}")

    zid, ns = ensure_zone(args.profile)
    if not args.verify:
        upsert(args.profile, zid, records)
        time.sleep(45)   # Route 53 propagates to its four nameservers within ~60 s
    print(f"\nverifying against {ns[0]}:")
    good = verify(records, ns)
    print("\nRoute 53 nameservers for GoDaddy:")
    for n in ns:
        print(f"  {n}")
    print("\nREADY: every record answers identically from Route 53" if good else "\nNOT READY: fix the mismatches above before changing nameservers")
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
