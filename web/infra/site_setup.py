#!/usr/bin/env python3
"""One-time, idempotent setup of the site's AWS resources (boto3; needs only a
profile with S3, CloudFront, ACM, Route 53 and STS permissions).

  certificate   ACM in us-east-1 (the only region CloudFront accepts) for
                openbiota.com and *.openbiota.com, DNS-validated through the
                Route 53 zone so it validates and renews on its own.
  bucket        private S3 bucket in us-west-1 (with the other site buckets),
                public access blocked, versioning on; CloudFront reads it
                through an Origin Access Control.
  distribution  CloudFront: HTTP/2+3, redirect to HTTPS, TLS 1.2 (2021),
                Brotli/gzip compression at the edge, security headers, 404
                page, PriceClass_All.
  function      the clean-URL CloudFront Function (infra/clean_urls.js:
                /index.html -> /, /docs/index.html -> /docs/, folder paths
                served from their index.html), created or updated to match
                the file, checked with CloudFront's test runner, published
                and attached to the distribution's viewer requests.
  dns           apex + www ALIAS A/AAAA records to the distribution.

    .venv/bin/python infra/site_setup.py --profile personal

Re-run until it reports every step complete; the certificate step waits for
validation, which takes a few minutes once the zone is live. Writes the ids
it creates to deploy.env for deploy_site.sh.

    .venv/bin/python infra/site_setup.py --profile personal --clean-urls

does the function step alone - deploy_site.sh runs it on every deploy, so an
edit to clean_urls.js goes live with the site and nothing else is touched.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import uuid
from pathlib import Path

import boto3
from botocore.exceptions import ClientError

DOMAIN = "openbiota.com"
BUCKET = "openbiota.com"
BUCKET_REGION = "us-west-1"
CERT_REGION = "us-east-1"
ZONE_ID = "Z00267273SHTHSY5M2IEY"
ENV_FILE = Path(__file__).resolve().parents[1] / "deploy.env"
CLOUDFRONT_HOSTED_ZONE = "Z2FDTNDATAQYW2"   # fixed id for every CloudFront distribution alias
#: AWS managed cache policy "CachingOptimized": honours origin Cache-Control, compression-aware keys
CACHING_OPTIMIZED = "658327ea-f89d-4fab-a63d-7e88639e58f6"
FUNCTION_NAME = f"{DOMAIN.replace(chr(46), chr(45))}-clean-urls"
FUNCTION_SOURCE = Path(__file__).with_name("clean_urls.js")
#: What the function must do with a request path (and query) before it is published:
#: a new path to serve, or the Location of a redirect.
FUNCTION_CHECKS = (
    ("/", {}, {"uri": "/index.html"}),
    ("/index.html", {}, {"location": "/"}),
    ("/index.html", {"utm_source": "hn"}, {"location": "/?utm_source=hn"}),
    ("/docs/index.html", {}, {"location": "/docs/"}),
    ("/docs", {}, {"location": "/docs/"}),
    ("/docs/", {}, {"uri": "/docs/index.html"}),
    ("/docs/install.html", {}, {"uri": "/docs/install.html"}),
    ("/research.html", {}, {"uri": "/research.html"}),
    ("/style.css", {"v": "202610071500"}, {"uri": "/style.css"}),
    ("/llms.txt", {}, {"uri": "/llms.txt"}),
)


class Clients:
    def __init__(self, profile: str) -> None:
        session = boto3.Session(profile_name=profile)
        self.s3 = session.client("s3", region_name=BUCKET_REGION)
        self.cf = session.client("cloudfront")
        self.acm = session.client("acm", region_name=CERT_REGION)
        self.r53 = session.client("route53")
        self.sts = session.client("sts")


def read_env() -> dict[str, str]:
    env: dict[str, str] = {}
    if ENV_FILE.is_file():
        for line in ENV_FILE.read_text().splitlines():
            if "=" in line and not line.startswith("#"):
                k, v = line.split("=", 1)
                env[k.strip()] = v.strip().strip('"')
    return env


def write_env(env: dict[str, str]) -> None:
    lines = ["# Deployment targets for deploy_site.sh. Nothing here is a secret.",
             "# Written by infra/site_setup.py; edit by hand to point at other resources, or to add",
             "# DOCS_SOURCE=\"../openbiota-gut-health-test\" (the software checkout, relative to this folder)",
             "# once the site has its own repository."]
    lines += [f'{k}="{v}"' for k, v in sorted(env.items())]
    ENV_FILE.write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
def certificate(c: Clients, env: dict[str, str]) -> str | None:
    arn = env.get("CERTIFICATE_ARN")
    if not arn:
        for cert in c.acm.list_certificates(MaxItems=100).get("CertificateSummaryList", []):
            if cert["DomainName"] == DOMAIN and f"*.{DOMAIN}" in (cert.get("SubjectAlternativeNameSummaries") or []):
                arn = cert["CertificateArn"]
        if not arn:
            arn = c.acm.request_certificate(
                DomainName=DOMAIN, SubjectAlternativeNames=[f"*.{DOMAIN}"], ValidationMethod="DNS",
                KeyAlgorithm="RSA_2048", IdempotencyToken="openbiotasite1",
                Tags=[{"Key": "site", "Value": DOMAIN}])["CertificateArn"]
            print(f"certificate requested: {arn}")
        env["CERTIFICATE_ARN"] = arn
        write_env(env)
    desc = c.acm.describe_certificate(CertificateArn=arn)["Certificate"]
    # the validation CNAME (shared by the apex and the wildcard) goes in the zone
    seen: set[str] = set()
    changes = []
    for opt in desc.get("DomainValidationOptions", []):
        rr = opt.get("ResourceRecord")
        if rr and rr["Name"] not in seen:
            seen.add(rr["Name"])
            changes.append({"Action": "UPSERT", "ResourceRecordSet": {
                "Name": rr["Name"], "Type": rr["Type"], "TTL": 300, "ResourceRecords": [{"Value": rr["Value"]}]}})
    if changes:
        c.r53.change_resource_record_sets(HostedZoneId=ZONE_ID, ChangeBatch={"Changes": changes})
        print(f"validation record in Route 53: {sorted(seen)}")
    print(f"certificate status: {desc['Status']}")
    return arn if desc["Status"] == "ISSUED" else None


def bucket(c: Clients) -> None:
    names = {b["Name"] for b in c.s3.list_buckets().get("Buckets", [])}
    if BUCKET not in names:
        c.s3.create_bucket(Bucket=BUCKET, CreateBucketConfiguration={"LocationConstraint": BUCKET_REGION})
        print(f"bucket created: s3://{BUCKET} ({BUCKET_REGION})")
    else:
        print(f"bucket exists: s3://{BUCKET}")
    c.s3.put_public_access_block(Bucket=BUCKET, PublicAccessBlockConfiguration={
        "BlockPublicAcls": True, "IgnorePublicAcls": True, "BlockPublicPolicy": True, "RestrictPublicBuckets": True})
    c.s3.put_bucket_versioning(Bucket=BUCKET, VersioningConfiguration={"Status": "Enabled"})
    # keep the last 10 versions of each object for 30 days: a bad deploy is recoverable, the bucket never grows without bound
    c.s3.put_bucket_lifecycle_configuration(Bucket=BUCKET, LifecycleConfiguration={"Rules": [{
        "ID": "trim-old-versions", "Status": "Enabled", "Filter": {"Prefix": ""},
        "NoncurrentVersionExpiration": {"NoncurrentDays": 30, "NewerNoncurrentVersions": 10}}]})
    c.s3.put_bucket_tagging(Bucket=BUCKET, Tagging={"TagSet": [{"Key": "site", "Value": DOMAIN}]})
    print("bucket: public access blocked, versioning on, old versions trimmed after 30 days")


def oac(c: Clients, env: dict[str, str]) -> str:
    if env.get("OAC_ID"):
        return env["OAC_ID"]
    for it in c.cf.list_origin_access_controls().get("OriginAccessControlList", {}).get("Items", []):
        if it["Name"] in (f"{DOMAIN}-s3", f"{DOMAIN.replace(chr(46), chr(45))}-s3"):
            env["OAC_ID"] = it["Id"]
            write_env(env)
            return it["Id"]
    res = c.cf.create_origin_access_control(OriginAccessControlConfig={
        "Name": f"{DOMAIN.replace(chr(46), chr(45))}-s3", "Description": f"CloudFront -> s3://{BUCKET}", "SigningProtocol": "sigv4",
        "SigningBehavior": "always", "OriginAccessControlOriginType": "s3"})
    env["OAC_ID"] = res["OriginAccessControl"]["Id"]
    write_env(env)
    print(f"origin access control created: {env['OAC_ID']}")
    return env["OAC_ID"]


def response_headers_policy(c: Clients) -> str:
    name = f"{DOMAIN.replace(chr(46), chr(45))}-security-headers"
    for p in c.cf.list_response_headers_policies(Type="custom").get("ResponseHeadersPolicyList", {}).get("Items", []):
        if p["ResponseHeadersPolicy"]["ResponseHeadersPolicyConfig"]["Name"] == name:
            return p["ResponseHeadersPolicy"]["Id"]
    res = c.cf.create_response_headers_policy(ResponseHeadersPolicyConfig={
        "Name": name, "Comment": "HSTS, nosniff, referrer policy, frame denial",
        "SecurityHeadersConfig": {
            "StrictTransportSecurity": {"Override": True, "IncludeSubdomains": True, "Preload": False, "AccessControlMaxAgeSec": 31536000},
            "ContentTypeOptions": {"Override": True},
            "FrameOptions": {"Override": True, "FrameOption": "DENY"},
            "ReferrerPolicy": {"Override": True, "ReferrerPolicy": "strict-origin-when-cross-origin"},
            "XSSProtection": {"Override": True, "Protection": True, "ModeBlock": True},
        }})
    print("response headers policy created")
    return res["ResponseHeadersPolicy"]["Id"]


def clean_urls_function(c: Clients) -> str:
    """The clean-URL function, brought level with infra/clean_urls.js: created if missing, updated if
    the code differs, run through FUNCTION_CHECKS on CloudFront's test runner, then published. Returns its ARN."""
    code = FUNCTION_SOURCE.read_bytes()
    config = {"Comment": "index.html out of the address bar; folder paths served from their index.html",
              "Runtime": "cloudfront-js-2.0"}
    try:
        desc = c.cf.describe_function(Name=FUNCTION_NAME, Stage="DEVELOPMENT")
        etag, arn = desc["ETag"], desc["FunctionSummary"]["FunctionMetadata"]["FunctionARN"]
        if c.cf.get_function(Name=FUNCTION_NAME, Stage="DEVELOPMENT")["FunctionCode"].read() != code:
            etag = c.cf.update_function(Name=FUNCTION_NAME, IfMatch=etag, FunctionConfig=config, FunctionCode=code)["ETag"]
            print("clean-url function updated from clean_urls.js")
        else:
            print("clean-url function matches clean_urls.js")
    except c.cf.exceptions.NoSuchFunctionExists:
        res = c.cf.create_function(Name=FUNCTION_NAME, FunctionConfig=config, FunctionCode=code)
        etag, arn = res["ETag"], res["FunctionSummary"]["FunctionMetadata"]["FunctionARN"]
        print(f"clean-url function created: {arn}")
    for uri, query, expected in FUNCTION_CHECKS:
        event = {"version": "1.0", "context": {"eventType": "viewer-request"}, "viewer": {"ip": "203.0.113.1"},
                 "request": {"method": "GET", "uri": uri, "querystring": {k: {"value": v} for k, v in query.items()},
                             "headers": {"host": {"value": DOMAIN}}, "cookies": {}}}
        result = c.cf.test_function(Name=FUNCTION_NAME, IfMatch=etag, Stage="DEVELOPMENT",
                                    EventObject=json.dumps(event).encode())["TestResult"]
        if result.get("FunctionErrorMessage"):
            sys.exit(f"clean-url function fails on {uri}: {result['FunctionErrorMessage']}")
        out = json.loads(result["FunctionOutput"])
        got = ({"location": out["response"]["headers"]["location"]["value"]} if out.get("response")
               else {"uri": out["request"]["uri"]})
        if got != expected:
            sys.exit(f"clean-url function: {uri} {query or ''} -> {got}, expected {expected}; not published")
    c.cf.publish_function(Name=FUNCTION_NAME, IfMatch=etag)
    print(f"clean-url function: {len(FUNCTION_CHECKS)} checks passed, published")
    return arn


def attach_function(c: Clients, dist_id: str, function_arn: str) -> None:
    """The clean-URL function on every viewer request of the distribution; nothing happens when it already is."""
    res = c.cf.get_distribution_config(Id=dist_id)
    config, etag = res["DistributionConfig"], res["ETag"]
    wanted = [{"FunctionARN": function_arn, "EventType": "viewer-request"}]
    if config["DefaultCacheBehavior"].get("FunctionAssociations", {}).get("Items") == wanted:
        print("distribution: clean-url function attached")
        return
    config["DefaultCacheBehavior"]["FunctionAssociations"] = {"Quantity": 1, "Items": wanted}
    c.cf.update_distribution(Id=dist_id, IfMatch=etag, DistributionConfig=config)
    print("distribution: clean-url function attached to viewer requests - deploying, ~5 min")


def distribution(c: Clients, env: dict[str, str], cert_arn: str, oac_id: str, function_arn: str) -> tuple[str, str]:
    if env.get("DISTRIBUTION_ID"):
        d = c.cf.get_distribution(Id=env["DISTRIBUTION_ID"])["Distribution"]
        print(f"distribution exists: {d['Id']} ({d['DomainName']}) status {d['Status']}")
        attach_function(c, d["Id"], function_arn)
        return d["Id"], d["DomainName"]
    rhp_id = response_headers_policy(c)
    config = {
        "CallerReference": str(uuid.uuid4()),
        "Comment": f"{DOMAIN} static site",
        "Enabled": True,
        "Aliases": {"Quantity": 2, "Items": [DOMAIN, f"www.{DOMAIN}"]},
        "DefaultRootObject": "index.html",
        "HttpVersion": "http2and3",
        "IsIPV6Enabled": True,
        "PriceClass": "PriceClass_All",
        "Origins": {"Quantity": 1, "Items": [{
            "Id": "s3-site", "DomainName": f"{BUCKET}.s3.{BUCKET_REGION}.amazonaws.com",
            "OriginAccessControlId": oac_id, "S3OriginConfig": {"OriginAccessIdentity": ""}}]},
        "DefaultCacheBehavior": {
            "TargetOriginId": "s3-site", "ViewerProtocolPolicy": "redirect-to-https",
            "AllowedMethods": {"Quantity": 2, "Items": ["GET", "HEAD"], "CachedMethods": {"Quantity": 2, "Items": ["GET", "HEAD"]}},
            "Compress": True, "CachePolicyId": CACHING_OPTIMIZED, "ResponseHeadersPolicyId": rhp_id,
            "FunctionAssociations": {"Quantity": 1, "Items": [{"FunctionARN": function_arn, "EventType": "viewer-request"}]},
        },
        "CustomErrorResponses": {"Quantity": 2, "Items": [
            {"ErrorCode": 403, "ResponsePagePath": "/404.html", "ResponseCode": "404", "ErrorCachingMinTTL": 10},
            {"ErrorCode": 404, "ResponsePagePath": "/404.html", "ResponseCode": "404", "ErrorCachingMinTTL": 10}]},
        "ViewerCertificate": {"ACMCertificateArn": cert_arn, "SSLSupportMethod": "sni-only", "MinimumProtocolVersion": "TLSv1.2_2021"},
    }
    d = c.cf.create_distribution(DistributionConfig=config)["Distribution"]
    env["DISTRIBUTION_ID"], env["DISTRIBUTION_DOMAIN"] = d["Id"], d["DomainName"]
    write_env(env)
    print(f"distribution created: {d['Id']} ({d['DomainName']}) - deploying, ~5 min")
    account = c.sts.get_caller_identity()["Account"]
    c.s3.put_bucket_policy(Bucket=BUCKET, Policy=json.dumps({"Version": "2012-10-17", "Statement": [{
        "Sid": "AllowCloudFrontServicePrincipalReadOnly", "Effect": "Allow",
        "Principal": {"Service": "cloudfront.amazonaws.com"}, "Action": "s3:GetObject",
        "Resource": f"arn:aws:s3:::{BUCKET}/*",
        "Condition": {"StringEquals": {"AWS:SourceArn": f"arn:aws:cloudfront::{account}:distribution/{d['Id']}"}}}]}))
    print("bucket policy: readable only by this distribution")
    return d["Id"], d["DomainName"]


def dns(c: Clients, cf_domain: str) -> None:
    changes = [{"Action": "UPSERT", "ResourceRecordSet": {
        "Name": f"{name}.", "Type": rtype,
        "AliasTarget": {"HostedZoneId": CLOUDFRONT_HOSTED_ZONE, "DNSName": f"{cf_domain}.", "EvaluateTargetHealth": False}}}
        for name in (DOMAIN, f"www.{DOMAIN}") for rtype in ("A", "AAAA")]
    c.r53.change_resource_record_sets(HostedZoneId=ZONE_ID, ChangeBatch={"Changes": changes})
    print(f"dns: {DOMAIN} and www -> {cf_domain} (A + AAAA alias)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default="personal")
    ap.add_argument("--clean-urls", action="store_true",
                    help="only publish the clean-URL function and attach it to the existing distribution (what deploy_site.sh runs)")
    args = ap.parse_args()
    c = Clients(args.profile)
    env = read_env()
    if args.clean_urls:
        if not env.get("DISTRIBUTION_ID"):
            sys.exit("no DISTRIBUTION_ID in deploy.env - run the full setup first")
        try:
            attach_function(c, env["DISTRIBUTION_ID"], clean_urls_function(c))
        except ClientError as exc:
            sys.exit(f"cloudfront: {exc}")
        return 0
    env.setdefault("AWS_PROFILE", args.profile)
    env.setdefault("BUCKET", BUCKET)
    env.setdefault("BUCKET_REGION", BUCKET_REGION)
    env.setdefault("ZONE_ID", ZONE_ID)
    write_env(env)

    bucket(c)
    oac_id = oac(c, env)
    cert_arn = certificate(c, env)
    if cert_arn is None:
        print("\nwaiting for the certificate to validate (its CNAME is in Route 53; the zone must be live at the registrar for ACM to see it)")
        for _ in range(40):
            time.sleep(15)
            status = c.acm.describe_certificate(CertificateArn=env["CERTIFICATE_ARN"])["Certificate"]["Status"]
            if status == "ISSUED":
                cert_arn = env["CERTIFICATE_ARN"]
                print("certificate ISSUED")
                break
            print(f"  … {status}")
        if cert_arn is None:
            print("certificate not issued yet; re-run after the nameservers point at Route 53")
            return 2
    try:
        function_arn = clean_urls_function(c)
        _dist_id, cf_domain = distribution(c, env, cert_arn, oac_id, function_arn)
    except ClientError as exc:
        sys.exit(f"cloudfront: {exc}")
    dns(c, cf_domain)
    print(f"\nsetup complete. deploy with: ./deploy_site.sh   (targets in {ENV_FILE})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
