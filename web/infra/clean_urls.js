// openbiota.com - clean URLs at the edge. A CloudFront Function on viewer
// requests (runtime cloudfront-js-2.0), published and attached to the
// distribution by infra/site_setup.py.
//
// The site is plain files in a private S3 bucket, which has no notion of a
// directory index below the root, and its pages link to one another by file
// name so that the folder also works from disk and on any host. This keeps
// "index.html" out of the address bar all the same:
//
//   /index.html         -> 301 /            the query string is kept; a #fragment
//   /docs/index.html    -> 301 /docs/       never reaches the edge and survives
//   /docs               -> 301 /docs/       the redirect in every browser
//   /                   -> served from /index.html       (no redirect)
//   /docs/              -> served from /docs/index.html
//
// Anything with a file extension - research.html, docs/install.html,
// style.css, llms.txt, favicon.ico - passes through untouched.

function handler(event) {
  var request = event.request;
  var uri = request.uri;
  if (uri.endsWith('/index.html')) {
    return redirect(uri.slice(0, -'index.html'.length) + queryString(request.querystring));
  }
  var name = uri.slice(uri.lastIndexOf('/') + 1);
  if (name !== '' && name.indexOf('.') === -1) {
    return redirect(uri + '/' + queryString(request.querystring));
  }
  if (uri.endsWith('/')) {
    request.uri = uri + 'index.html';
  }
  return request;
}

function redirect(location) {
  return {
    statusCode: 301,
    statusDescription: 'Moved Permanently',
    headers: {
      location: { value: location },
      'cache-control': { value: 'public, max-age=86400' }
    }
  };
}

// The query string as it arrived, "?a=1&b=2" or "" - parameters are passed
// through as received, repeated ones included, without decoding or re-encoding.
function queryString(params) {
  var parts = [];
  for (var name in params) {
    var values = params[name].multiValue || [params[name]];
    for (var i = 0; i < values.length; i++) {
      parts.push(values[i].value === '' ? name : name + '=' + values[i].value);
    }
  }
  return parts.length ? '?' + parts.join('&') : '';
}
