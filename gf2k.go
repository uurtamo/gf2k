package "github.com/uurtamo/gf2k"

import (

	"math/bits"
	"fmt"
)

type Field struct {

	// for now, these are all base 2

	// Size will be k for GF(2^k)
	// in particular, this means that poly will require 1 extra bit
	// of storage above the value of k

	Size uint64          
	Poly *Polynomial
	Generator *Element
}

type Polynomial struct {

	// for now, these will all be primitive trinomials

	Packed []uint64
	HighWord uint64
	HighBit uint64
	LowWord uint64
	LowBit uint64

	IsPrimitive bool
}

type Element struct {
	
	Packed []uint64
	Field *Field
}

func (E *Element) Init (F *Field) {

	E.Field  = F
	E.Packed = make([]uint64, (F.Size/64)+1)
}

func (F *Field) Init (High uint64, Low uint64, Primitive bool, generator *Element) {

	F.Size=High

	poly := new(Polynomial)

	highbit := High
	lowbit := Low
	highword := highbit/64
	lowword := lowbit/64
	
	poly.HighWord=highword
	poly.HighBit=highbit
	poly.LowWord=lowword
	poly.LowBit=lowbit
	poly.IsPrimitive=Primitive

	// fmt.Printf("packbits size is %d\n",highword+1)

	packbits := make([]uint64, highword+1)

	packbits[0]=1
	packbits[lowbit/64]|=1<<(poly.LowBit % 64)
	packbits[highbit/64]|=1<<(poly.HighBit % 64)

	if Primitive==true {
		poly.IsPrimitive=true

		prim:=new(Element)
		prim.Init(F)
		prim.Packed[0]=2
		F.Generator=prim
	} else {
		poly.IsPrimitive=false
		F.Generator=nil
	}

	poly.Packed=packbits
	F.Poly=poly

	if Primitive==false {F.Generator=generator}
}

func (e1 *Element) Copy(e2 *Element) {

	//	for i:=0; i<=e2.Field.Poly.HighWord; i++ {
	//		e1.packed[i]=e2.packed[i]
	//	}
	copy(e1.Packed,e2.Packed)
}

func Multiply(e1 *Element, e2 *Element) *Element {

	e3:=new(Element)
	e3.Init(e1.Field)

	e1Copy:=new(Element)
	e1Copy.Init(e1.Field)

	e2Copy:=new(Element)
	e2Copy.Init(e1.Field)

	e1Copy.Copy(e1)
	e2Copy.Copy(e2)

	PeasantAlgorithm(e1Copy,e2Copy,e3,e1.Field)
	return e3
}

func MultiplySparseFaster(e1 *Element, e2 *Element) *Element {

	var setBits uint64
	var wordCount uint64
	var iter uint64
	var LL []uint64

	e3:=new(Element)
	e3.Init(e1.Field)

	for wordCount=0; wordCount<=e2.Field.Poly.HighWord; wordCount++ {

		if bits.OnesCount64(e2.Packed[wordCount])>0 {

			setBits+=uint64(bits.OnesCount64(e2.Packed[wordCount]))

			for iter=0; iter<=63; iter++ {
				if ((1 << iter) & e2.Packed[wordCount]) > 0 {
					LL=append(LL,(wordCount*64)+iter)
				}
			}
		}
	}

	fmt.Printf("setbits in second argument is %d\n",setBits)

	fmt.Printf("bits are at %v\n",LL)

	NaiveMultiplyFaster(e1,LL,e3,e1.Field)
	return e3
}

func PeasantAlgorithm(e1 *Element, e2 *Element, e3 *Element, F *Field) {

	var i uint64
	var j uint64
	var bit uint64
	var carry uint64

	// e1 * e2 will be stored in e3

	for i=0; i<=F.Poly.HighBit; i++ {

		// if the bottom bit of e2 is a 1, XOR the product by e1

		if (e2.Packed[0] & 1) == 1 {
			for j=0; j<= F.Poly.HighWord; j++ {
				e3.Packed[j] ^= e1.Packed[j]
			}
		}

		//		if bits.TrailingZeros64(e2.Packed[0]) == 0 {
		//			for j=0; j<= F.Poly.HighWord; j++ {
		//				e3.Packed[j] ^= e1.Packed[j]
		//			}
		//		}

		// shift e2 right by 1 bit

		for j=0; j< F.Poly.HighWord; j++ {
			e2.Packed[j] >>= 1
			bit = (e2.Packed[j+1] & 1)
			e2.Packed[j] |= (bit <<63)
		}

		e2.Packed[F.Poly.HighWord] >>= 1;

		// shift e1 left by 1 bit

		for j=F.Poly.HighWord; j>0; j-- {
			e1.Packed[j] <<= 1
			bit = ((e1.Packed[j-1] & (1<<63)) >> 63)
			e1.Packed[j] |= bit
		}

		e1.Packed[0] <<= 1

		// carry is the top (valid) bit of e1
		// if (carry == 1), XOR e1 by the field poly
		// this just means that the representation has overflown


		carry = e1.Packed[F.Poly.HighWord] & (1 << (F.Poly.HighBit % 64))

		if (carry > 0) {

			e1.Packed[F.Poly.HighWord] ^= F.Poly.Packed[F.Poly.HighWord]

			if F.Poly.LowWord != F.Poly.HighWord {
				e1.Packed[F.Poly.LowWord] ^= F.Poly.Packed[F.Poly.LowWord]
			}

			if (F.Poly.LowWord > 0) {
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

	fmt.Printf("multiplying %v by bits in %v for field %v\n",e1,e2,F)

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

		WordOffset   = i/64

		if i % 64 == 0 {

			// we really just need to slide and XOR

			for j=0; j<=F.Poly.HighWord; j++ {
				doublewide[ j + WordOffset] ^= e1.Packed[ j ]
			}
		} else {
			// we need to half-slide and XOR two part-words at a time
			// we'll need a low mask and a high mask:

			lowmask  = (1<<(64-i))-1           // the 64-i lowest bits of a word
			highmask = ((1<<i)-1) << (64-i)   // the i(th) highest bits of a word

			// say we have [w3] [w2] [w1] [w0] ... and we need
			// to shift left by 17 bits. we'll end up with
			// [47 bits of 0's, 17 high bits of w3]
			// [47 low bits of w3, 17 high bits of w2]
			// [47 low bits of w2, 17 high bits of w1]
			// [47 low bits of w1, 17 high bits of w0]
			// [47 low bits of w0, 17 bits of 0's]

			// first the highest target word

			dump = 0
			dump |= ((highmask & e1.Packed[F.Poly.HighWord]) >> (64-i))
			doublewide[ F.Poly.HighWord + WordOffset ] ^= dump

			// now the lowest target word

			dump = 0
			dump |= ((lowmask & e1.Packed[0]) << i)
			doublewide[ WordOffset ] ^= dump

			// there is a case where we might have a very small field where
			// we can't do this loop, so we need to check for it

			if F.Poly.HighWord > 0 {
				for j=1; j<=(F.Poly.HighWord-1); j++ {

					dump     = 0
					dump |= ((lowmask & e1.Packed[j]) << i)
					dump |= ((highmask & e1.Packed[j-1]) >> (64-i))
					doublewide[ j + WordOffset] ^= dump
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

		fmt.Printf("bitptr: %d F.Poly.HighBit %d\n",bitptr,F.Poly.HighBit);

		// everything to the left of bitptr should already be
		// a zero

		trail = (bitptr % 64) + 1
		thunderbird = bitptr/64
		window=0

		// trail is the leftmost portion of the window that we should obtain
		// from the current word expressed as a bitlength. if it has a leading one,
		// we'll zero it out along with the ones for the rest of the field polynomial, then
		// we'll count the leading zeros, and advance bitptr that many
		// slots

		if (trail == 64) || (bitptr < 64) {

			window = doublewide[thunderbird] << (64-trail)

		} else {

			window |= ((((1<<trail)-1) & doublewide[thunderbird]) << (64-trail))
			window |= (((((1<<(64-trail))-1) << trail) & doublewide[thunderbird-1]) >> trail)
		}

		skippo=bits.LeadingZeros64(window)

		if skippo==0 {

			// mod out by the field polynomial

			doublewide[thunderbird] ^= (1<<(trail -1))
			doublewide[(bitptr-gap1)/64] ^= (1 << ((bitptr-gap1) % 64))
			doublewide[(bitptr-gap1-gap2)/64] ^= (1 << ((bitptr-gap1-gap2) % 64))

		} else {

			bitptr -= uint64(skippo)
		}
	}

	// either bitptr is exactly at F.Poly.HighBit or it's past it. if it's
	// at it, then we need to check to see if that's a 1 and if so mod by
	// the field polynomial.

	if bitptr == F.Poly.HighBit {

		if ((1<<(bitptr % 64)) & doublewide[bitptr/64]) > 0 {
			doublewide[bitptr/64] ^= (1 << (bitptr % 64))
			doublewide[(bitptr-gap1)/64] ^= (1 << ((bitptr-gap1) % 64))
			doublewide[(bitptr-gap1-gap2)/64] ^= (1 << ((bitptr-gap1-gap2) % 64))
		}
	}

	// now we need to pack the remainder into the return object
	
	for i=0; i<=F.Poly.HighWord; i++ {

		e3.Packed[i]=doublewide[i]
	}
}



func Add(e1 *Element, e2 *Element) *Element {

	var i uint64

	e3:=new(Element)
	e3.Init(e1.Field)

	hw:=e3.Field.Poly.HighWord

	for i=0; i<=hw; i++ {
		e3.Packed[i]=e1.Packed[i]^e2.Packed[i]
	}
	return(e3)
}

func Invert(e1 *Element) *Element {

	// raise to the power (2^k)-2 OR
	// extended euclidean algorithm
	// should be a linear number of squares and multiplies (in k) in the first case

	var i uint64

	result:=new(Element)
	result.Init(e1.Field)
	result.Packed[0]=1

	e2:=new(Element)
	e2.Init(e1.Field)
	e2.Copy(e1)

	hw:=e1.Field.Poly.HighBit

	// 2^k    == 100000..000
	// 2^k -1 == 011111..111
	// 2^k -2 == 011111..110


	for i=1; i<=(hw-1); i++ {
		//		fmt.Printf("e2 before sq\n")
		//		e2.Print()
		e2=Square(e2)
		//		fmt.Printf("e2 after sq\n")
		//		e2.Print()
		result=Multiply(e2,result)
		//		fmt.Printf("result after multiply\n")
		//		result.Print()
	}
	e2=result
	return(e2)
}

func Square(e1 *Element) *Element {

	e3:=new(Element)
	e3.Init(e1.Field)

	e3=Multiply(e1,e1)
	return(e3)
}

func (e1 *Element) MultiplyInPlace(e2 *Element) {

	e3:=new(Element)
	e3.Init(e1.Field)

	e2Copy:=new(Element)
	e2Copy.Init(e1.Field)

	e2Copy.Copy(e2)

	PeasantAlgorithm(e1,e2Copy,e3,e1.Field)

	//	fmt.Printf(" e1 before\n")
	//	e1.Print()

	e1=e3

	//	fmt.Printf(" e1 after\n")
	//	e1.Print()
}

func (e1 *Element) SquareInPlace() {

	e3:=new(Element)
	e3.Init(e1.Field)

	e2Copy:=new(Element)
	e2Copy.Init(e1.Field)

	e2Copy.Copy(e1)

	PeasantAlgorithm(e1,e2Copy,e3,e1.Field)
	e1=e3
}

func (e1 *Element) AddInPlace(e2 *Element) {

	var i uint64

	for i=0; i<=e1.Field.Poly.HighWord; i++ {
		e1.Packed[i]^=e2.Packed[i]
	}
}

func (e1 *Element) InvertInPlace() {

}

func (e1 *Element) HighBitIsSet() bool{

	if (e1.Packed[e1.Field.Poly.HighWord] >> (e1.Field.Poly.HighBit % 64))==1 {
		return true
	} else {
		return false
	}
}

func (e1 *Element) Overflow() bool{

	if (e1.Packed[e1.Field.Poly.HighWord] >> (e1.Field.Poly.HighBit % 64))>1 {
		return true
	} else {
		return false
	}
}

func (e1 *Element) Print() {

	var i uint64

	for i=0; i<=e1.Field.Poly.HighWord; i++ {

		fmt.Printf("%v ",e1.Packed[i])
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

	hw:=e1.Field.Poly.HighWord
	lw:=e1.Field.Poly.LowWord

	e1.Packed[hw]^=e1.Field.Poly.Packed[hw]
	e1.Packed[lw]^=e1.Field.Poly.Packed[lw]

	if (lw>0) {
		e1.Packed[0]^=e1.Field.Poly.Packed[0]
	}
}

