import pytest

from sunpy.util.multimethod import MultiMethod, FAIL

def test_super():
    class String(str):
        pass
    
    mm = MultiMethod(lambda *a: a)
    
    @mm.add_dec(str, str)
    def foo(foo, bar):
        return 'String'
    
    @mm.add_dec(String, str)
    def foo(foo, bar):
        return 'Fancy', mm.super(super(String, foo), bar)
    
    assert mm('foo', 'bar') == 'String'
    assert mm(String('foo'), 'bar') == ('Fancy', 'String')


def test_override():
    class String(str):
        pass
    
    mm = MultiMethod(lambda *a: a)
    
    @mm.add_dec(str, str)
    def foo(foo, bar):
        return 'String'
    
    pytest.raises(
        TypeError, mm.add_dec(String, str, override=FAIL), lambda x, y: None
    )