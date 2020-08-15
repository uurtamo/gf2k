package gf2k

import (
	"errors"
	"fmt"
	"math/bits"
	"os"
)

type Field struct {

	// for now, these are all base 2

	// Size will be k for GF(2^k)
	// in particular, this means that poly will require 1 extra bit
	// of storage above the value of k

	Size      uint64
	Poly      *Polynomial
	Generator *Element
}

type Polynomial struct {

	// for now, these will all be primitive trinomials

	Packed   []uint64
	HighWord uint64
	HighBit  uint64
	LowWord  uint64
	LowBit   uint64

	IsPrimitive bool
}

type Element struct {
	Packed []uint64
	Field  *Field
}

func (E *Element) Init(F *Field) {

	E.Field = F
	E.Packed = make([]uint64, (F.Size/64)+1)
}

func (F *Field) Init(High uint64, Low uint64, Primitive bool, generator *Element) {

	F.Size = High

	poly := new(Polynomial)

	highbit := High
	lowbit := Low
	highword := highbit / 64
	lowword := lowbit / 64

	poly.HighWord = highword
	poly.HighBit = highbit
	poly.LowWord = lowword
	poly.LowBit = lowbit
	poly.IsPrimitive = Primitive

	// fmt.Printf("packbits size is %d\n",highword+1)

	packbits := make([]uint64, highword+1)

	packbits[0] = 1
	packbits[lowbit/64] |= 1 << (poly.LowBit % 64)
	packbits[highbit/64] |= 1 << (poly.HighBit % 64)

	if Primitive == true {
		poly.IsPrimitive = true

		prim := new(Element)
		prim.Init(F)
		prim.Packed[0] = 2
		F.Generator = prim
	} else {
		poly.IsPrimitive = false
		F.Generator = nil
	}

	poly.Packed = packbits
	F.Poly = poly

	if Primitive == false {
		F.Generator = generator
	}
}

func IsZero(e1 *Element) bool {

	var cn uint64
	var bitsCn uint64

	cn = 0

	for bitsCn == 0 && cn <= e1.Field.Poly.HighWord {

		bitsCn += uint64(bits.OnesCount64(e1.Packed[cn]))
		cn++
	}

	if bitsCn > 0 {
		return false
	} else {
		return true
	}
}

func (e1 *Element) Copy(e2 interface{}) {

	//	for i:=0; i<=e2.Field.Poly.HighWord; i++ {
	//		e1.packed[i]=e2.packed[i]
	//	}

	switch e2 := e2.(type) {

	// if e2 is a [] uint64

	case []uint64:

		copy(e1.Packed, e2)

	// if e2 is a * Element

	case *Element:

		copy(e1.Packed, e2.Packed)
	}
}

func Multiply(e1 *Element, e2 *Element) *Element {

	e3 := new(Element)
	e3.Init(e1.Field)

	e1Copy := new(Element)
	e1Copy.Init(e1.Field)

	e2Copy := new(Element)
	e2Copy.Init(e1.Field)

	e1Copy.Copy(e1)
	e2Copy.Copy(e2)

	PeasantAlgorithm(e1Copy, e2Copy, e3, e1.Field)
	return e3
}

func MultiplySparse(e1 *Element, e2 *Element) *Element {

	var setBits uint64
	var wordCount uint64
	var iter uint64
	var LL []uint64

	e3 := new(Element)
	e3.Init(e1.Field)

	for wordCount = 0; wordCount <= e2.Field.Poly.HighWord; wordCount++ {

		if bits.OnesCount64(e2.Packed[wordCount]) > 0 {

			setBits += uint64(bits.OnesCount64(e2.Packed[wordCount]))

			for iter = 0; iter <= 63; iter++ {
				if ((1 << iter) & e2.Packed[wordCount]) > 0 {
					LL = append(LL, (wordCount*64)+iter)
				}
			}
		}
	}

	//	fmt.Printf("setbits in second argument is %d\n",setBits)

	//	fmt.Printf("bits are at %v\n",LL)

	NaiveMultiplyFaster(e1, LL, e3, e1.Field)
	return e3
}

func PeasantAlgorithm(e1 *Element, e2 *Element, e3 *Element, F *Field) {

	var i uint64
	var j uint64
	var bit uint64
	var carry uint64

	// e1 * e2 will be stored in e3

	for i = 0; i <= F.Poly.HighBit; i++ {

		// if the bottom bit of e2 is a 1, XOR the product by e1

		if (e2.Packed[0] & 1) == 1 {
			for j = 0; j <= F.Poly.HighWord; j++ {
				e3.Packed[j] ^= e1.Packed[j]
			}
		}

		//		if bits.TrailingZeros64(e2.Packed[0]) == 0 {
		//			for j=0; j<= F.Poly.HighWord; j++ {
		//				e3.Packed[j] ^= e1.Packed[j]
		//			}
		//		}

		// shift e2 right by 1 bit

		for j = 0; j < F.Poly.HighWord; j++ {
			e2.Packed[j] >>= 1
			bit = (e2.Packed[j+1] & 1)
			e2.Packed[j] |= (bit << 63)
		}

		e2.Packed[F.Poly.HighWord] >>= 1

		// shift e1 left by 1 bit

		for j = F.Poly.HighWord; j > 0; j-- {
			e1.Packed[j] <<= 1
			bit = ((e1.Packed[j-1] & (1 << 63)) >> 63)
			e1.Packed[j] |= bit
		}

		e1.Packed[0] <<= 1

		// carry is the top (valid) bit of e1
		// if (carry == 1), XOR e1 by the field poly
		// this just means that the representation has overflown

		carry = e1.Packed[F.Poly.HighWord] & (1 << (F.Poly.HighBit % 64))

		if carry > 0 {

			e1.Packed[F.Poly.HighWord] ^= F.Poly.Packed[F.Poly.HighWord]

			if F.Poly.LowWord != F.Poly.HighWord {
				e1.Packed[F.Poly.LowWord] ^= F.Poly.Packed[F.Poly.LowWord]
			}

			if F.Poly.LowWord > 0 {
				e1.Packed[0] ^= F.Poly.Packed[0]
			}
		}
	}

	//	fmt.Printf("multiply complete\n")

	//	fmt.Printf(" e1\n")
	//	e1.Print()

	//	fmt.Printf(" e2\n")
	//	e2.Print()

	//	fmt.Printf(" result:\n")
	//	e3.Print()
}

func NaiveMultiplyFaster(e1 *Element, e2 []uint64, e3 *Element, F *Field) {

	//	fmt.Printf("multiplying %v by bits in %v for field %v\n",e1,e2,F)

	// in the special case where only a constant number of bits are
	// going to be set in the second argument, we can just pass the
	// offsets of those bits and do the multiplication in linear time

	// we'll do this in 2x the space of the first argument, naively
	// multiplying, then reduce the answer with the field polynomial
	// in one pass, storing the result in the third argument as with
	// the peasant algorithm

	// we'll rotate through the product through adjacent zeroes, using
	// the leading zeroes builtin function, (which executes in ~4 clocks
	// on haswell) to skip unneccesary runs of computation. this should
	// save roughly

	// this means occasionally checking the number of leading zeroes
	// of upcoming words

	// we could do this with multibyte lookup tables, but the builtin
	// leading zeroes function is fast enough to make that an unnecessary
	// complexity

	var i uint64
	var j uint64
	var highmask uint64
	var lowmask uint64
	var WordOffset uint64
	var dump uint64
	var bitptr uint64
	var window uint64
	var thunderbird uint64
	var gap1 uint64
	var gap2 uint64
	var trail uint64
	var skippo int

	// eventually we'll be pulling objects like this out of a pool in order
	// to shave a few operations off of the allocation

	doublewide := make([]uint64, 2*(F.Poly.HighWord+1))

	// e1.Packed looks like: [b63 b62...b2 b1 b0] [b63 b62...b2 b1 b0] .. [b63 b62...b2 b1 b0]
	// e2 looks like: [o1] [o2] .. [oh] for offsets 1 through h

	for _, i = range e2 {

		WordOffset = i / 64

		if i%64 == 0 {

			// we really just need to slide and XOR

			for j = 0; j <= F.Poly.HighWord; j++ {
				doublewide[j+WordOffset] ^= e1.Packed[j]
			}
		} else {
			// we need to half-slide and XOR two part-words at a time
			// we'll need a low mask and a high mask:

			lowmask = (1 << (64 - i)) - 1         // the 64-i lowest bits of a word
			highmask = ((1 << i) - 1) << (64 - i) // the i(th) highest bits of a word

			// say we have [w3] [w2] [w1] [w0] ... and we need
			// to shift left by 17 bits. we'll end up with
			// [47 bits of 0's, 17 high bits of w3]
			// [47 low bits of w3, 17 high bits of w2]
			// [47 low bits of w2, 17 high bits of w1]
			// [47 low bits of w1, 17 high bits of w0]
			// [47 low bits of w0, 17 bits of 0's]

			// first the highest target word

			dump = 0
			dump |= ((highmask & e1.Packed[F.Poly.HighWord]) >> (64 - i))
			doublewide[F.Poly.HighWord+WordOffset] ^= dump

			// now the lowest target word

			dump = 0
			dump |= ((lowmask & e1.Packed[0]) << i)
			doublewide[WordOffset] ^= dump

			// there is a case where we might have a very small field where
			// we can't do this loop, so we need to check for it

			if F.Poly.HighWord > 0 {
				for j = 1; j <= (F.Poly.HighWord - 1); j++ {

					dump = 0
					dump |= ((lowmask & e1.Packed[j]) << i)
					dump |= ((highmask & e1.Packed[j-1]) >> (64 - i))
					doublewide[j+WordOffset] ^= dump
				}
			}
		}
	}

	// doublewide now contains the product, unreduced, that we'd like
	// to reduce into singlewide form

	// we can skip forward over any leading 0's at any time as long
	// as we have at least enough bits left to fill up a singlewide

	bitptr = 2 * F.Poly.HighBit

	gap1 = F.Poly.HighBit - F.Poly.LowBit
	gap2 = F.Poly.LowBit

	// we're going to construct a virtual window of 64 bits starting
	// at bitptr and check that window for leading zeros, then skip
	// over them

	for bitptr > F.Poly.HighBit {

		//		fmt.Printf("bitptr: %d F.Poly.HighBit %d\n",bitptr,F.Poly.HighBit);

		// everything to the left of bitptr should already be
		// a zero

		trail = (bitptr % 64) + 1
		thunderbird = bitptr / 64
		window = 0

		// trail is the leftmost portion of the window that we should obtain
		// from the current word expressed as a bitlength. if it has a leading one,
		// we'll zero it out along with the ones for the rest of the field polynomial, then
		// we'll count the leading zeros, and advance bitptr that many
		// slots

		if (trail == 64) || (bitptr < 64) {

			window = doublewide[thunderbird] << (64 - trail)

		} else {

			window |= ((((1 << trail) - 1) & doublewide[thunderbird]) << (64 - trail))
			window |= (((((1 << (64 - trail)) - 1) << trail) & doublewide[thunderbird-1]) >> trail)
		}

		skippo = bits.LeadingZeros64(window)

		if skippo == 0 {

			// mod out by the field polynomial

			doublewide[thunderbird] ^= (1 << (trail - 1))
			doublewide[(bitptr-gap1)/64] ^= (1 << ((bitptr - gap1) % 64))
			doublewide[(bitptr-gap1-gap2)/64] ^= (1 << ((bitptr - gap1 - gap2) % 64))

		} else {

			if uint64(skippo) > bitptr {
				bitptr = 0

			} else {
				bitptr -= uint64(skippo)
			}
		}
	}

	// either bitptr is exactly at F.Poly.HighBit or it's past it. if it's
	// at it, then we need to check to see if that's a 1 and if so mod by
	// the field polynomial.

	if bitptr == F.Poly.HighBit {

		if ((1 << (bitptr % 64)) & doublewide[bitptr/64]) > 0 {
			doublewide[bitptr/64] ^= (1 << (bitptr % 64))
			doublewide[(bitptr-gap1)/64] ^= (1 << ((bitptr - gap1) % 64))
			doublewide[(bitptr-gap1-gap2)/64] ^= (1 << ((bitptr - gap1 - gap2) % 64))
		}
	}

	// now we need to pack the remainder into the return object

	for i = 0; i <= F.Poly.HighWord; i++ {

		e3.Packed[i] = doublewide[i]
	}
}

func Add(e1 *Element, e2 *Element) *Element {

	var i uint64

	e3 := new(Element)
	e3.Init(e1.Field)

	hw := e3.Field.Poly.HighWord

	for i = 0; i <= hw; i++ {
		e3.Packed[i] = e1.Packed[i] ^ e2.Packed[i]
	}
	return (e3)
}

func Invert(e1 *Element) *Element {

	// raise to the power (2^k)-2 OR
	// extended euclidean algorithm
	// should be a linear number of squares and multiplies (in k) in the first case

	//	var i uint64

	//	result:=new(Element)
	//	result.Init(e1.Field)
	//	result.Packed[0]=1

	e2 := new(Element)
	e2.Init(e1.Field)
	e2.ForceReduce()

	e3 := new(Element)
	e3.Init(e1.Field)

	//	hw:=e1.Field.Poly.HighBit

	// 2^k    == 100000..000
	// 2^k -1 == 011111..111
	// 2^k -2 == 011111..110

	//	for i=1; i<=(hw-1); i++ {
	//		fmt.Printf("e2 before sq\n")
	//		e2.Print()
	//		e2=Square(e2)
	//		fmt.Printf("e2 after sq\n")
	//		e2.Print()
	//		result=MultiplySparse(e2,result)
	//		fmt.Printf("result after multiply\n")
	//		result.Print()
	//	}
	//	e2=result

	e3 = ExtendedEuclideanAlgorithm(e2, e1)

	return (e3)
}

func Square(e1 *Element) *Element {

	e3 := new(Element)
	e3.Init(e1.Field)

	e3 = MultiplySparse(e1, e1)
	return (e3)
}

func (e1 *Element) MultiplyInPlace(e2 *Element) {

	e3 := new(Element)
	e3.Init(e1.Field)

	e2Copy := new(Element)
	e2Copy.Init(e1.Field)

	e2Copy.Copy(e2)

	PeasantAlgorithm(e1, e2Copy, e3, e1.Field)

	//	fmt.Printf(" e1 before\n")
	//	e1.Print()

	e1 = e3

	//	fmt.Printf(" e1 after\n")
	//	e1.Print()
}

func (e1 *Element) SquareInPlace() {

	e3 := new(Element)
	e3.Init(e1.Field)

	e2Copy := new(Element)
	e2Copy.Init(e1.Field)

	e2Copy.Copy(e1)

	PeasantAlgorithm(e1, e2Copy, e3, e1.Field)
	e1 = e3
}

func (e1 *Element) AddInPlace(e2 *Element) {

	var i uint64

	for i = 0; i <= e1.Field.Poly.HighWord; i++ {
		e1.Packed[i] ^= e2.Packed[i]
	}
}

func (e1 *Element) InvertInPlace() {

	e2 := new(Element)
	e2.Init(e1.Field)

	e2.ForceReduce()

	e1 = ExtendedEuclideanAlgorithm(e2, e1)

}

func (e1 *Element) HighBitIsSet() bool {

	if (e1.Packed[e1.Field.Poly.HighWord] >> (e1.Field.Poly.HighBit % 64)) == 1 {
		return true
	} else {
		return false
	}
}

func (e1 *Element) Overflow() bool {

	if (e1.Packed[e1.Field.Poly.HighWord] >> (e1.Field.Poly.HighBit % 64)) > 1 {
		return true
	} else {
		return false
	}
}

func (e1 *Element) Print() {

	var i uint64

	for i = 0; i <= e1.Field.Poly.HighWord; i++ {

		fmt.Printf("%v ", e1.Packed[i])
	}
	fmt.Printf("\n")
}

func (e1 *Element) Reduce() {

	if e1.HighBitIsSet() {
		e1.ForceReduce()
	}
}

func (e1 *Element) ForceReduce() {

	//	for (i:=0; i<=e1.Field.Poly.HighWord; i++) {
	//		e1[i]^=e1.Field.Poly.packed[i]
	//	}

	hw := e1.Field.Poly.HighWord
	lw := e1.Field.Poly.LowWord

	e1.Packed[hw] ^= e1.Field.Poly.Packed[hw]
	e1.Packed[lw] ^= e1.Field.Poly.Packed[lw]

	if lw > 0 {
		e1.Packed[0] ^= e1.Field.Poly.Packed[0]
	}
}

func FindHighBit(array []uint64, start uint64) (uint64, error) {

	var wordmax uint64
	var bitptr uint64
	var maxBit uint64
	var highOne uint64
	var trail uint64
	var thunderbird uint64
	var window uint64
	var skippo uint64
	var skipped uint64

	// words are offsets into the array, so indexed starting at zero

	wordmax = uint64(len(array) - 1)

	// bits are counted starting at zero

	bitptr = ((wordmax + 1) * 64) - 1
	maxBit = bitptr

	if (array[wordmax] & (1 << 63)) > 0 {

		highOne = maxBit

	} else {

		if start > 0 {

			maxBit = start
			bitptr = start
		}

		highOne = bitptr

		skipped = 0

		for highOne >= maxBit {

			// trail is indexed starting at 1
			// bitptr is indexed starting at 0

			trail = (bitptr % 64) + 1
			thunderbird = bitptr / 64

			if (trail == 64) || (bitptr < 64) {

				window = array[thunderbird] << (64 - trail)

			} else {

				window |= ((((1 << trail) - 1) & array[thunderbird]) << (64 - trail))
				window |= array[thunderbird-1] >> trail
			}

			skippo = uint64(bits.LeadingZeros64(window))

			skipped += skippo

			if skippo == 0 {
				highOne = bitptr

			} else {

				if skipped > maxBit {
					return 0, errors.New("FindHighBit: Was passed the all zeros array\n")

				} else {
					bitptr -= uint64(skippo)
				}
			}
		}
	}
	return highOne, nil
}

func ModularReduction(modulus []uint64, item []uint64) []uint64 {

	//	fmt.Printf("ModularReduction: modulus: %v item: %v\n",modulus,item)

	quotient := make([]uint64, len(item))
	remainder := make([]uint64, len(modulus))

	var highBitItem uint64
	var highBitModulus uint64

	var err error

	highBitItem, err = FindHighBit(item, 0)

	// the item is zero

	if err != nil {
		return item
	}

	highBitModulus, err = FindHighBit(modulus, 0)

	// the modulus is zero

	if err != nil {
		return modulus
	}

	if len(item) < len(modulus) {

		//		fmt.Printf("ModularReduction: remainder: %v\n",item)
		return item

	} else if highBitItem < highBitModulus {

		// the item is already itself mod the modulus

		//		fmt.Printf("ModularReduction: remainder: %v\n",item)
		return item

	} else {

		//		fmt.Printf("ModularReduction calling PolynomialDivide: divisor: %v dividend: %v quotient %v remainder %v\n",modulus, item, quotient, remainder)

		PolynomialDivide(modulus, item, quotient, remainder)

		//		fmt.Printf("ModularReduction: remainder: %v\n",remainder)
		return (remainder)
	}
}

func PolynomialDivide(divisor []uint64, dividend []uint64, quotient []uint64, remainder []uint64) {

	// we are going to divide the dividend by the divisor, storing the result in
	// the quotient and the remainder

	// we assume that quotient and remainder are passed in already zeroed

	//	fmt.Printf("PolynomialDivide: input: divisor %v dividend %v quotient %v remainder %v \n",divisor,dividend,quotient,remainder)

	var bitptr uint64
	var trail uint64
	var thunderbird uint64
	var window uint64
	var skippo uint64
	var skipped uint64

	var divisorHighBit uint64
	var divisorHighWord uint64

	var dividendHighBit uint64
	var dividendHighWord uint64

	var offsetGap uint64

	var wordCounter uint64

	var divisorWindow uint64

	var dividendOffset uint64
	var divisorOffset uint64

	var wordGap uint64
	var rem uint64

	var err error

	// first we need to determine the bitwidth of the dividend and the divisor

	divisorHighBit, err = FindHighBit(divisor, 0)
	if err != nil {

		fmt.Printf("PolynomialDivide: attempt to divide by 0\n")
		os.Exit(13)
	}

	dividendHighBit, err = FindHighBit(dividend, 0)
	if err != nil {
		fmt.Printf("PolynomialDivide: dividend is 0; this is likely not what you intended\n")
		os.Exit(13)
	}

	//	fmt.Printf("PolynomialDivide: divisorHighBit %v, dividendHighBit %v\n",divisorHighBit, dividendHighBit)

	divisorOffset = divisorHighBit % 64
	dividendOffset = dividendHighBit % 64

	if divisorHighBit > dividendHighBit {

		// we need to reduce the divisor by the dividend,
		// meaning take the remainder after dividing the
		// divisor by the dividend

		divisor = ModularReduction(dividend, divisor)
		divisorHighBit, err = FindHighBit(divisor, divisorHighBit)
		if err != nil {
			fmt.Printf("PolynomialDivide after calling ModularReduction: divisor is 0. this is likely not what you intended\n")
			os.Exit(13)
		}
		divisorOffset = divisorHighBit % 64
	}

	bitptr = dividendHighBit

	//	fmt.Printf("PolynomialDivide: bitptr %v, divisorHighBit %v\n",bitptr, divisorHighBit)

	rem = 1

	for (bitptr >= divisorHighBit) && (rem != 0) {

		// everything to the left of bitptr should already be
		// a zero

		trail = (bitptr % 64) + 1
		thunderbird = bitptr / 64
		window = 0

		// trail is the leftmost portion of the window that we should obtain
		// from the current word expressed as a bitlength. if it has a leading one,
		// we'll zero it out along with the ones for the rest of the field polynomial, then
		// we'll count the leading zeros, and advance bitptr that many
		// slots

		if (trail == 64) || (bitptr < 64) {

			window = dividend[thunderbird] << (64 - trail)
			//			fmt.Printf("PolynomialDivide: trail == 64 (%v) || bitptr < 64 (%v), window %v\n",trail, bitptr, window)

		} else {

			window |= ((((1 << trail) - 1) & dividend[thunderbird]) << (64 - trail))
			window |= (((((1 << (64 - trail)) - 1) << trail) & dividend[thunderbird-1]) >> trail)
			//			fmt.Printf("PolynomialDivide: trail %v bitptr %v window %v\n",trail, bitptr, window)
		}

		skippo = uint64(bits.LeadingZeros64(window))

		//		fmt.Printf("PolynomialDivide: skippo %v\n",skippo)

		skipped += skippo

		if skippo == 0 {

			// we've found the first set 1 bit in the leftmost nonzero word of the dividend

			// now we need to mod out by the divisor -- a 1 will be stored in the quotient
			// at this offset

			divisorWindow = 0

			dividendHighWord = bitptr / 64
			divisorHighWord = divisorHighBit / 64

			dividendOffset = bitptr % 64

			quotient[dividendHighWord] |= (1 << (dividendOffset))

			wordGap = (dividendHighWord) - (divisorHighWord)

			//			fmt.Printf("PolynomialDivide: dividendOffset: %v, dividendHighWord: %v, quotient: %v, wordGap: %v\n",dividendOffset, dividendHighWord, quotient, wordGap)

			if dividendOffset == divisorOffset {

				//				fmt.Printf("PolynomialDivide: dividendOffset %v, divisorOffset %v\n",dividendOffset,divisorOffset);

				// the divisor and the dividend have the same alignment, so XORing
				// doesn't require any shifting

				for wordCounter = wordGap; wordCounter <= dividendHighWord; wordCounter++ {
					dividend[wordCounter] ^= divisor[wordCounter-wordGap]
					//					fmt.Printf("PolynomialDivide: wordCounter %v, dividend[wc] %v\n",wordCounter,dividend[wordCounter])
				}

				//				fmt.Printf("PolynomialDivide: dividendOffset==divisorOffset, dividend %v\n",dividend)

			} else if dividendOffset > divisorOffset {

				//				fmt.Printf("PolynomialDivide: dividendOffset %v, divisorOffset %v\n",dividendOffset,divisorOffset);

				offsetGap = dividendOffset - divisorOffset

				for wordCounter = dividendHighWord; wordCounter >= wordGap+1; wordCounter-- {

					// this is the case where there are at least two words required
					// to describe the divisor

					// either the dividend is aligned, or neither are aligned.
					// all we know so far is that they have different
					// alignment and that the dividend is 'more left' in their respective
					// first nonzero words

					divisorWindow |= divisor[wordCounter-wordGap] << offsetGap
					divisorWindow |= ((((1 << offsetGap) - 1) << (64 - offsetGap)) & divisor[wordCounter-1]) >> (64 - offsetGap)

					dividend[wordCounter] ^= divisorWindow

					//					fmt.Printf("PolynomialDivide: wordCounter %v, divisorWindow %v, dividend[wc] %v\n",wordCounter,divisorWindow,dividend[wordCounter])
				}

				// now we handle the bottommost word of the divisor

				dividend[wordGap] ^= (divisor[0] << offsetGap)
				//				fmt.Printf("PolynomialDivide: wordGap %v, divisor[0] %v, offsetGap %v\n",wordGap,divisor[0],offsetGap)

			} else {

				//				fmt.Printf("PolynomialDivide: dividendOffset %v, divisorOffset %v\n",dividendOffset,divisorOffset);

				// dividend offset should be < divsor offset

				// first we handle the topmost word of the divisor

				offsetGap = divisorOffset - dividendOffset

				dividend[dividendHighWord] ^= (divisor[divisorHighWord] >> offsetGap)

				// note -- divisorHighBit should never be bigger than dividendHighBit
				// when *both* are 1 word long

				if (divisorHighWord == 0 && dividendHighWord == 0) || (wordGap == 0) {

					fmt.Printf("this should never happen\n")
					fmt.Printf("divisorOffset: %v dividendOffset: %v divisorHighWord: %v wordGap: %v\n", divisorOffset, dividendOffset, divisorHighWord, wordGap)
					fmt.Printf("divisorHighBit: %v dividendHighBit: %v\n", divisorHighBit, dividendHighBit)
					os.Exit(13)
				}

				for wordCounter = dividendHighWord - 1; wordCounter >= wordGap; wordCounter-- {

					// this is the case where there are at least two words required
					// to describe the divisor

					// either the divisor is aligned, or neither are aligned.
					// all we know so far is that they have different
					// alignment and that the dividend is 'more right' in their respective
					// first nonzero words

					divisorWindow |= divisor[wordCounter-wordGap] >> offsetGap
					divisorWindow |= (((1 << offsetGap) - 1) & divisor[wordCounter+1]) << (64 - offsetGap)

					dividend[wordCounter] ^= divisorWindow
				}
			}

		} else {

			if skippo > bitptr {

				// this will happen if dividend is all zeros, for instance, because leading
				// zeros will be 64
				bitptr = 0
				rem = 0

			} else {

				bitptr -= uint64(skippo)
			}
		}
	}
	for wordCounter = 0; wordCounter <= divisorHighWord; wordCounter++ {

		remainder[wordCounter] = dividend[wordCounter]
	}

	//	fmt.Printf("PolynomialDivide quotient %v\n",quotient)
	//	fmt.Printf("PolynomialDivide remainder %v\n",remainder)
}

func ExtendedEuclideanAlgorithm(e1 *Element, e2 *Element) *Element {

	var quotient []uint64
	var remainder []uint64

	Pass := new(Element)
	Pass.Init(e1.Field)

	Pass2 := new(Element)
	Pass2.Init(e1.Field)

	R := new(Element)
	R.Init(e1.Field)

	NewR := new(Element)
	NewR.Init(e1.Field)

	S := new(Element)
	S.Init(e1.Field)

	NewS := new(Element)
	NewS.Init(e1.Field)

	T := new(Element)
	T.Init(e1.Field)

	NewT := new(Element)
	NewT.Init(e1.Field)

	S.Packed[0] = 1
	NewS.Packed[0] = 0

	T.Packed[0] = 0
	NewT.Packed[0] = 1

	// to invert e2 in field e1:

	R.Packed[e1.Field.Poly.HighWord] ^= e1.Field.Poly.Packed[e1.Field.Poly.HighWord]

	if e1.Field.Poly.LowWord != e1.Field.Poly.HighWord {
		R.Packed[e1.Field.Poly.LowWord] ^= e1.Field.Poly.Packed[e1.Field.Poly.LowWord]
	}
	if e1.Field.Poly.LowWord > 0 {
		R.Packed[0] ^= e1.Field.Poly.Packed[0]
	}

	NewR.Copy(e2)

	//		fmt.Printf("ExtendedEuclideanAlgorithm\n")
	//		fmt.Printf("NewR %v\n",NewR)
	//		fmt.Printf("R %v\n",R)

	for !IsZero(NewR) {

		//				fmt.Printf("Calling PolynomialDivide for NewR %v and R %v\n",NewR.Packed,R.Packed)

		quotient = make([]uint64, e1.Field.Poly.HighWord+1)
		remainder = make([]uint64, e1.Field.Poly.HighWord+1)

		PolynomialDivide(NewR.Packed, R.Packed, quotient, remainder)

		R.Copy(NewR)
		NewR.Copy(remainder)

		Pass2.Copy(NewT)

		Pass.Copy(quotient)
		NewT = Add(T, Multiply(Pass, NewT))

		NewT.Reduce()

		T.Copy(Pass2)

		//				fmt.Printf("   Now we have R %v, NewR %v, T %v, and NewT %v\n",R.Packed,NewR.Packed,T.Packed,NewT.Packed)

	}

	// T should contain the inverse of e2

	//	fmt.Printf("   T %v, should be the inverse \n",T.Packed)

	return (T)
}
