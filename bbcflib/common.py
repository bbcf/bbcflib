###############################################################################
def normalize_url(url):
    """Produce a fixed form for an HTTP URL.

    Make sure the URL begins with http:// and does *not* end with a /.

    >>> normalize_url('http://www.google.com')
    'http://www.google.com'
    >>> normalize_url('http://www.google.com/')
    'http://www.google.com'
    >>> normalize_url('www.google.com/')
    'http://www.google.com'
    """
    url = url.lower()
    if not(url.startswith("http://")):
        url = "http://" + url
    if url.endswith("/"):
        url = url[:-1]
    return url

###############################################################################
def natural_sort(item):
    '''
    Will sort strings that contain numbers correctly

    >>> l = ['v1.3.12', 'v1.3.3', 'v1.2.5', 'v1.2.15', 'v1.2.3', 'v1.2.1']
    >>> l.sort(key=natural_sort)
    >>> l.__repr__()
    "['v1.2.1', 'v1.2.3', 'v1.2.5', 'v1.2.15', 'v1.3.3', 'v1.3.12']"
    '''
    import re
    def try_int(s):
        try: return int(s)
        except: return s
    return map(try_int, re.findall(r'(\d+|\D+)', item))
