# Copyright (c) 2011 Florian Mayer <florian.mayer@bitsrc.org>

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from warnings import warn
from itertools import izip

SILENT = 0
WARN = 1
FAIL = 2

def _fmt_t(types):
    return ', '.join(type_.__name__ for type_ in types)

class MultiMethod(object):
    def __init__(self, get):
        self.get = get
        
        self.methods = []
        self.cache = {}
    
    def add(self, fun, types, override=SILENT):
        overriden = False
        if override:
            for signature, _ in self.methods:
                if all(issubclass(a, b) for a, b in izip(types, signature)):
                    overriden = True
        if overriden and override == FAIL:
            raise TypeError
        elif overriden and override == WARN:
            warn(
                'Definition (%s) overrides prior definition (%s).' %
                (_fmt_t(types), _fmt_t(signature)),
                stacklevel=3
            )
        elif overriden:
            raise ValueError('Invalid value for override.')
        self.methods.append((types, fun))
    
    def add_dec(self, *types, **kwargs):
        self.cache = {}
        def _dec(fun):
            self.add(fun, types, kwargs.get('override', SILENT))
            return fun
        return _dec
    
    def __call__(self, *args, **kwargs):
        objs = self.get(*args, **kwargs)
        
        types = tuple(map(type, objs))
        
        # This code is duplicate for performace reasons.
        cached = self.cache.get(types, None)
        if cached is not None:
            return cached(*args, **kwargs)
        
        for signature, fun in reversed(self.methods):
            if all(issubclass(ty, sig) for ty, sig in zip(types, signature)):
                self.cache[types] = fun
                return fun(*args, **kwargs)
        raise TypeError
    
    def super(self, *args, **kwargs):
        objs = self.get(*args, **kwargs)
        types = tuple(
            [
                x.__thisclass__.__mro__[1] if isinstance(x, super) else type(x)
                for x in objs
            ]
        )
        nargs = [
            x.__self__ if isinstance(x, super) else x
            for x in args
        ]
        
        for k, elem in kwargs.iteritems():
            if isinstance(elem, super):
                kwargs[k] = elem.__self__
        
        # This code is duplicate for performace reasons.
        cached = self.cache.get(types, None)
        if cached is not None:
            return cached(*nargs, **kwargs)
        
        for signature, fun in reversed(self.methods):
            if all(issubclass(ty, sig) for ty, sig in zip(types, signature)):
                self.cache[types] = fun
                return fun(*nargs, **kwargs)
        raise TypeError


if __name__ == '__main__':
    class String(str):
        def foo(self, foo, bar):
            return 'String', foo, bar
    
    
    mm = MultiMethod(lambda *a: a)
    
    @mm.add_dec(str, str)
    def foo(foo, bar):
        return 'String'
    
    @mm.add_dec(String, str, override=WARN)
    def foo(foo, bar):
        return 'Fancy', mm.super(super(String, foo), bar)
        
    @mm.add_dec(int, str)
    def foo(foo, bar):
        return 'Int - String'
    
    assert mm('foo', 'bar') == 'String'
    assert mm(1, 'bar') == 'Int - String'
    
    @mm.add_dec(int, int)
    def foo(foo, bar):
        return foo + bar
    
    assert mm(1, 2) == 3
    
    assert mm(String('foo'), 'bar') == ('Fancy', 'String')
    from time import time
    
    s = time()
    for _ in xrange(int(1e6)):
        mm('foo', 'bar')
    print time() - s

    s = time()
    for _ in xrange(int(1e6)):
        foo('foo', 'bar')
    print time() - s

    st = String('foo')
    s = time()
    for _ in xrange(int(1e6)):
        st.foo('foo', 'bar')
    print time() - s


    s = time()
    for _ in xrange(int(1e6)):
        mm.super(super(String, String('foo')), 'bar')
    print time() - s
