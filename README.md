
Implementacja kodu BCH
======================

###Wymagania###

- python 2.7
- biblioteka numpy


##Uwaga##

Aby umożliwić stworzenie tablicy dodawania dla ciała GF(2^8) należy zmodyfikować ciało funkcji polydiv z biblioteki numpy.polynomial:

```python

[c1, c2] = pu.as_series([c1, c2])
if c2[-1] == 0 :
    raise ZeroDivisionError()

len1 = len(c1)
len2 = len(c2)
if len2 == 1 :
    return c1/c2[-1], c1[:1]*0
elif len1 < len2 :
    return c1[:1]*0, c1
else :
    dlen = len1 - len2
    scl = c2[-1]
    c2  = c2[:-1]/scl
    i = dlen
    j = len1 - 1
    while i >= 0 :
        c1[i:j] -= (c2*c1[j])%2  # dodane % 2
        c1[i:j] = c1[i:j]%2  # dodane % 2
        i -= 1
        j -= 1
    return c1[j+1:]/scl, pu.trimseq(c1[:j+1])
```
