#!/bin/bash
# Deploy openbiota.com/ to S3 + CloudFront.
#
#   ./deploy_site.sh            build, confirm, publish the edge function, upload, invalidate
#   ./deploy_site.sh --dry-run  build and show what would be uploaded
#   ./deploy_site.sh --no-docs  skip the documentation build (landing-page edits only)
#
# Targets come from deploy.env (written by infra/site_setup.py; not secret).
# The documentation Markdown comes from the software repository's checkout,
# DOCS_SOURCE (environment, or the DOCS_SOURCE line of deploy.env, relative to
# this folder: "../openbiota-gut-health-test" for a sibling checkout); see README.md.
#
# Cache policy: HTML, XML, TXT and the manifest are `no-cache` so browsers
# revalidate them on every visit and pick up new asset versions; CSS, JS,
# fonts and images are cached for a year because their URLs carry a version
# stamp (CSS/JS) or change when the file does (fonts, images are replaced,
# not edited). CloudFront compresses text at the edge (Brotli/gzip), so
# nothing is pre-compressed here.
set -euo pipefail

WEB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SITE="$WEB/openbiota.com"
ENV_FILE="$WEB/deploy.env"
PY="${PYTHON:-$WEB/.venv/bin/python}"
[ -x "$PY" ] || PY="python3"

DRY_RUN=""
NO_DOCS=""
for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN="--dryrun" ;;
    --no-docs) NO_DOCS="--no-docs" ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

[ -f "$ENV_FILE" ] || { echo "no $ENV_FILE - run infra/site_setup.py first" >&2; exit 1; }
# shellcheck disable=SC1090
source "$ENV_FILE"
: "${AWS_PROFILE:?}" "${BUCKET:?}" "${BUCKET_REGION:?}"
if [ -z "$DRY_RUN" ]; then
  : "${DISTRIBUTION_ID:?missing in deploy.env - infra/site_setup.py has not finished (certificate still validating?)}"
fi

command -v aws >/dev/null || { echo "aws CLI not found (brew install awscli)" >&2; exit 1; }

echo "== build"
"$PY" "$WEB/build_site.py" $NO_DOCS

echo "== check"
# nothing personal, no secrets, no build leftovers
if grep -rIl -E "AKIA[0-9A-Z]{16}|BEGIN (RSA|OPENSSH) PRIVATE KEY" "$SITE" >/dev/null; then
  echo "refusing to deploy: a credential-like string is in the site" >&2; exit 1
fi
if grep -rIl -E "ERICK_|MONICA_|KRY799|FCG748|BGE623|BXP835|AMD614|FXX745" "$SITE" >/dev/null; then
  echo "refusing to deploy: a local sample identifier is in the site" >&2; exit 1
fi
find "$SITE" -name .DS_Store -delete

if [ -z "$DRY_RUN" ]; then
  echo
  echo "This uploads $SITE to s3://$BUCKET and invalidates CloudFront $DISTRIBUTION_ID (https://openbiota.com)."
  read -r -p "Deploy? [y/N] " answer
  [ "$answer" = "y" ] || { echo "not deployed"; exit 1; }
fi

if [ -z "$DRY_RUN" ]; then
  echo "== edge"
  # the clean-URL function (infra/clean_urls.js) is part of the site: published and
  # attached to the distribution here, so an edit to it ships with the next deploy;
  # a no-op when nothing changed
  "$PY" "$WEB/infra/site_setup.py" --profile "$AWS_PROFILE" --clean-urls
fi

echo "== stamp"
# a new URL is a new cache entry: browsers holding an old page fetch new assets immediately
STAMP="$(date +%Y%m%d%H%M)"
for page in "$SITE"/index.html "$SITE"/research.html; do
  [ -f "$page" ] && sed -i '' -E "s/(\.(css|js)\?v=)[0-9]+/\1${STAMP}/g" "$page"
done

SYNC=(aws s3 sync "$SITE/" "s3://$BUCKET/" --profile "$AWS_PROFILE" --region "$BUCKET_REGION" --delete
      --exclude "README.md" --exclude ".DS_Store" --exclude "*/.DS_Store" $DRY_RUN)

echo "== upload: long-lived assets"
"${SYNC[@]}" --exclude "*.html" --exclude "*.html.md" --exclude "*.xml" --exclude "*.txt" --exclude "*.webmanifest" --exclude "*.json" \
  --cache-control "public, max-age=31536000, immutable"

echo "== upload: revalidated documents"
"${SYNC[@]}" --exclude "*" --include "*.html" --include "*.html.md" --include "*.xml" --include "*.txt" --include "*.webmanifest" --include "*.json" \
  --cache-control "no-cache"

# content types the CLI does not infer
if [ -z "$DRY_RUN" ]; then
  aws s3 cp "s3://$BUCKET/site.webmanifest" "s3://$BUCKET/site.webmanifest" --profile "$AWS_PROFILE" --region "$BUCKET_REGION" \
    --content-type "application/manifest+json" --cache-control "no-cache" --metadata-directive REPLACE --only-show-errors
  for md in $(cd "$SITE" && find . -name "*.html.md" -o -name "llms.txt" | sed 's#^\./##'); do
    aws s3 cp "s3://$BUCKET/$md" "s3://$BUCKET/$md" --profile "$AWS_PROFILE" --region "$BUCKET_REGION" \
      --content-type "text/markdown; charset=utf-8" --cache-control "no-cache" --metadata-directive REPLACE --only-show-errors
  done
fi

if [ -z "$DRY_RUN" ]; then
  echo "== invalidate"
  aws cloudfront create-invalidation --distribution-id "$DISTRIBUTION_ID" --paths "/*" --profile "$AWS_PROFILE" \
    --query "Invalidation.{id:Id,status:Status}" --output text
  echo
  echo "deployed: https://openbiota.com  (CloudFront serves the new files within a minute or two)"
fi
