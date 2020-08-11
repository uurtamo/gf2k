# gf2k

Implements finite field arithmetic over GF(2^k)

This gives you the basic operations:

Multiply, Add, Invert

It is mostly useful for when k > 64, since then you're doing multiword arithmetic.

For very small values of k, lookup tables are faster.

For very large values of k, goroutines are useful but require specialized locking or lock-free setups to avoid word updates from clobbering one another.

The easiest but somewhat cumbersome lock-free method is to put all of the XORs into CAS loops. Slightly nicer is to set aside ranges of words using a lock-free range locking structure.

Some niceties to the single-threaded implementation:

* If one of the two arguments that you're multiplying is expected to have a small number of set bits, there's a specialized multiply routine that you can call which will reduce it to time linear in the number of set bits in that argument.

* If you are very short on space, there is a specialized multiply routine that uses space the size of its arguments.

Note that on Haswell and after, counting the number of leading or trailing zeros in a word takes a trivial amount of time, so is almost always worth doing if it means that another efficiency can be had as a result of it being greater than 1.

Similarly, on Haswell and after, counting the total number of bits set in a word is very efficient.

Efficiencies to come and to be tested: (TODO)

* lock-free range locks for word range locking
* CAS loops for XORs

Timings for various values of k: (TODO)
