from .mocks import MOCK_HASH


def test_cache_basic(cache):
    cache.download('http://example.com/abc.text')
    assert cache._downloader.times_called == 1


def test_cache_caching(cache):
    cache.download('http://example.com/abc.text')
    cache.download('http://example.com/abc.text')
    assert cache._downloader.times_called == 1


def test_cache_redownload(cache):
    cache.download('http://example.com/abc.text')
    cache.download('http://example.com/abc.text', redownload=True)
    assert cache._downloader.times_called == 2


def test_get_by_url(cache):
    cache.download('http://example.com/file_name')
    details = cache._get_by_url('http://example.com/file_name')
    assert details['file_path'].endswith('file_name')


def test_get_by_url_fail(cache):
    cache.download('http://example.com/file_name')
    details = cache._get_by_url('http://example.com/file')
    assert details is None


def test_get_by_hash(cache):
    cache.download('http://example.com/file_name')
    details = cache.get_by_hash(MOCK_HASH)
    assert details['file_path'].endswith('file_name')


def test_get_by_hash_fail(cache):
    cache.download('http://example.com/file_name')
    details = cache.get_by_hash('wrong_hash')
    assert details is None
