
#include <cstddef>
#include <iostream>
#include <string>

#include "../include/cadet/StringUtil.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Reference implementation https://github.com/veorq/SipHash

#include <stdio.h>
#include <string.h>

/* default: SipHash-2-4 */
#define MAXLEN 64
#define KEYLEN 16
#ifndef DOUBLE
#define HASHLEN 8
#else
#define HASHLEN 16
#endif

#define cROUNDS 2
#define dROUNDS 4

#define ROTL(x,b) (uint64_t)( ((x) << (b)) | ( (x) >> (64 - (b))) )

#define U32TO8_LE(p, v)                                         \
  (p)[0] = (uint8_t)((v)      ); (p)[1] = (uint8_t)((v) >>  8); \
  (p)[2] = (uint8_t)((v) >> 16); (p)[3] = (uint8_t)((v) >> 24);

#define U64TO8_LE(p, v)                        \
  U32TO8_LE((p),     (uint32_t)((v)      ));   \
  U32TO8_LE((p) + 4, (uint32_t)((v) >> 32));

#define U8TO64_LE(p)            \
  (((uint64_t)((p)[0])      ) | \
   ((uint64_t)((p)[1]) <<  8) | \
   ((uint64_t)((p)[2]) << 16) | \
   ((uint64_t)((p)[3]) << 24) | \
   ((uint64_t)((p)[4]) << 32) | \
   ((uint64_t)((p)[5]) << 40) | \
   ((uint64_t)((p)[6]) << 48) | \
   ((uint64_t)((p)[7]) << 56))

#define SIPROUND                                        \
  do {                                                  \
    v0 += v1; v1=ROTL(v1,13); v1 ^= v0; v0=ROTL(v0,32); \
    v2 += v3; v3=ROTL(v3,16); v3 ^= v2;                 \
    v0 += v3; v3=ROTL(v3,21); v3 ^= v0;                 \
    v2 += v1; v1=ROTL(v1,17); v1 ^= v2; v2=ROTL(v2,32); \
  } while(0)


#define TRACE

int  siphashRefMain( uint8_t *out, const uint8_t *in, uint64_t inlen, const uint8_t *k )
{
  /* "somepseudorandomlygeneratedbytes" */
  uint64_t v0 = 0x736f6d6570736575ULL;
  uint64_t v1 = 0x646f72616e646f6dULL;
  uint64_t v2 = 0x6c7967656e657261ULL;
  uint64_t v3 = 0x7465646279746573ULL;
  uint64_t b;
  uint64_t k0 = U8TO64_LE( k );
  uint64_t k1 = U8TO64_LE( k + 8 );
  uint64_t m;
  int i;
  const uint8_t *end = in + inlen - ( inlen % sizeof( uint64_t ) );
  const int left = inlen & 7;
  b = ( ( uint64_t )inlen ) << 56;
  v3 ^= k1;
  v2 ^= k0;
  v1 ^= k1;
  v0 ^= k0;

#ifdef DOUBLE
  v1 ^= 0xee;
#endif

  for ( ; in != end; in += 8 )
  {
    m = U8TO64_LE( in );
    v3 ^= m;

    TRACE;
    for( i=0; i<cROUNDS; ++i ) SIPROUND;

    v0 ^= m;
  }

  switch( left )
  {
  case 7: b |= ( ( uint64_t )in[ 6] )  << 48;
  case 6: b |= ( ( uint64_t )in[ 5] )  << 40;
  case 5: b |= ( ( uint64_t )in[ 4] )  << 32;
  case 4: b |= ( ( uint64_t )in[ 3] )  << 24;
  case 3: b |= ( ( uint64_t )in[ 2] )  << 16;
  case 2: b |= ( ( uint64_t )in[ 1] )  <<  8;
  case 1: b |= ( ( uint64_t )in[ 0] ); break;
  case 0: break;
  }


  v3 ^= b;

  TRACE;
  for( i=0; i<cROUNDS; ++i ) SIPROUND;

  v0 ^= b;

#ifndef DOUBLE
  v2 ^= 0xff;
#else
  v2 ^= 0xee;
#endif

  TRACE;
  for( i=0; i<dROUNDS; ++i ) SIPROUND;

  b = v0 ^ v1 ^ v2  ^ v3;
  U64TO8_LE( out, b );

#ifdef DOUBLE
  v1 ^= 0xdd;

  TRACE;
  for( i=0; i<dROUNDS; ++i ) SIPROUND;

  b = v0 ^ v1 ^ v2  ^ v3;
  U64TO8_LE( out+8, b );
#endif

  return 0;
}

inline uint64_t sipHashRef(char const* data, unsigned int len)
{
    uint8_t in[MAXLEN], out[HASHLEN], k[KEYLEN];

    for (int i = 0; i < KEYLEN; ++i)
        k[i] = 0;
    
    siphashRefMain( out, reinterpret_cast<const uint8_t*>(data), len, k );

    return *(reinterpret_cast<uint64_t*>(out));
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <typename Enumeration>
inline typename std::underlying_type<Enumeration>::type as_integer(Enumeration const value)
{
    return static_cast<typename std::underlying_type<Enumeration>::type>(value);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using cadet::util::SipHash24;
using cadet::util::SipHash24runtime;
using cadet::operator "" _hash;

// Force compile time evaluation
enum class Bla : uint64_t
{
	Hallo = "Hallo"_hash,
	DasIstEinLangerTestAuchMitZahlenDrin12 = "DasIstEinLangerTestAuchMitZahlenDrin12"_hash,
	A = "A"_hash,
	AB = "AB"_hash,
	ABC = "ABC"_hash,
	ABCD = "ABCD"_hash,
	ABCDE = "ABCDE"_hash,
	ABCDEF = "ABCDEF"_hash,
	ABCDEFG = "ABCDEFG"_hash,
	ABCDEFGH = "ABCDEFGH"_hash,
	ABCDEFGHI = "ABCDEFGHI"_hash,
	ABCDEFGHIJ = "ABCDEFGHIJ"_hash,
	ABCDEFGHIJK = "ABCDEFGHIJK"_hash,
	ABCDEFGHIJKL = "ABCDEFGHIJKL"_hash,
	ABCDEFGHIJKLM = "ABCDEFGHIJKLM"_hash,
	ABCDEFGHIJKLMN = "ABCDEFGHIJKLMN"_hash,
	ABCDEFGHIJKLMNO = "ABCDEFGHIJKLMNO"_hash,
	ABCDEFGHIJKLMNOP = "ABCDEFGHIJKLMNOP"_hash,
	ABCDEFGHIJKLMNOPQ = "ABCDEFGHIJKLMNOPQ"_hash,
	ABCDEFGHIJKLMNOPQR = "ABCDEFGHIJKLMNOPQR"_hash,
	ABCDEFGHIJKLMNOPQRS = "ABCDEFGHIJKLMNOPQRS"_hash,
	ABCDEFGHIJKLMNOPQRST = "ABCDEFGHIJKLMNOPQRST"_hash,
	ABCDEFGHIJKLMNOPQRSTU = "ABCDEFGHIJKLMNOPQRSTU"_hash,
	ABCDEFGHIJKLMNOPQRSTUV = "ABCDEFGHIJKLMNOPQRSTUV"_hash,
	ABCDEFGHIJKLMNOPQRSTUVW = "ABCDEFGHIJKLMNOPQRSTUVW"_hash,
	ABCDEFGHIJKLMNOPQRSTUVWX = "ABCDEFGHIJKLMNOPQRSTUVWX"_hash,
	ABCDEFGHIJKLMNOPQRSTUVWXY  = "ABCDEFGHIJKLMNOPQRSTUVWXY"_hash,
	ABCDEFGHIJKLMNOPQRSTUVWXYZ = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"_hash,
};


int main(int argc, char** argv)
{
	std::cout << "Hallo: " << sipHashRef("Hallo", 5) << std::endl;
	std::cout << "Hallo: " << SipHash24runtime("Hallo") << std::endl;
	std::cout << "Hallo: " << SipHash24("Hallo") << std::endl;
	std::cout << "Hallo: " << as_integer(Bla::Hallo) << std::endl;

	std::cout << "DasIstEinLangerTestAuchMitZahlenDrin12: " << sipHashRef("DasIstEinLangerTestAuchMitZahlenDrin12", 38) << std::endl;
	std::cout << "DasIstEinLangerTestAuchMitZahlenDrin12: " << SipHash24runtime("DasIstEinLangerTestAuchMitZahlenDrin12") << std::endl;
	std::cout << "DasIstEinLangerTestAuchMitZahlenDrin12: " << SipHash24("DasIstEinLangerTestAuchMitZahlenDrin12") << std::endl;
	std::cout << "DasIstEinLangerTestAuchMitZahlenDrin12: " << as_integer(Bla::DasIstEinLangerTestAuchMitZahlenDrin12) << std::endl;

	std::cout << "A: " << sipHashRef("A", 1) << std::endl;
	std::cout << "A: " << SipHash24runtime("A") << std::endl;
	std::cout << "A: " << SipHash24("A") << std::endl;
	std::cout << "A: " << as_integer(Bla::A) << std::endl;

	std::cout << "AB: " << sipHashRef("AB", 2) << std::endl;
	std::cout << "AB: " << SipHash24runtime("AB") << std::endl;
	std::cout << "AB: " << SipHash24("AB") << std::endl;
	std::cout << "AB: " << as_integer(Bla::AB) << std::endl;

	std::cout << "ABC: " << sipHashRef("ABC", 3) << std::endl;
	std::cout << "ABC: " << SipHash24runtime("ABC") << std::endl;
	std::cout << "ABC: " << SipHash24("ABC") << std::endl;
	std::cout << "ABC: " << as_integer(Bla::ABC) << std::endl;

	std::cout << "ABCD: " << sipHashRef("ABCD", 4) << std::endl;
	std::cout << "ABCD: " << SipHash24runtime("ABCD") << std::endl;
	std::cout << "ABCD: " << SipHash24("ABCD") << std::endl;
	std::cout << "ABCD: " << as_integer(Bla::ABCD) << std::endl;

	std::cout << "ABCDE: " << sipHashRef("ABCDE", 5) << std::endl;
	std::cout << "ABCDE: " << SipHash24runtime("ABCDE") << std::endl;
	std::cout << "ABCDE: " << SipHash24("ABCDE") << std::endl;
	std::cout << "ABCDE: " << as_integer(Bla::ABCDE) << std::endl;

	std::cout << "ABCDEF: " << sipHashRef("ABCDEF", 6) << std::endl;
	std::cout << "ABCDEF: " << SipHash24runtime("ABCDEF") << std::endl;
	std::cout << "ABCDEF: " << SipHash24("ABCDEF") << std::endl;
	std::cout << "ABCDEF: " << as_integer(Bla::ABCDEF) << std::endl;

	std::cout << "ABCDEFG: " << sipHashRef("ABCDEFG", 7) << std::endl;
	std::cout << "ABCDEFG: " << SipHash24runtime("ABCDEFG") << std::endl;
	std::cout << "ABCDEFG: " << SipHash24("ABCDEFG") << std::endl;
	std::cout << "ABCDEFG: " << as_integer(Bla::ABCDEFG) << std::endl;

	std::cout << "ABCDEFGH: " << sipHashRef("ABCDEFGH", 8) << std::endl;
	std::cout << "ABCDEFGH: " << SipHash24runtime("ABCDEFGH") << std::endl;
	std::cout << "ABCDEFGH: " << SipHash24("ABCDEFGH") << std::endl;
	std::cout << "ABCDEFGH: " << as_integer(Bla::ABCDEFGH) << std::endl;

	std::cout << "ABCDEFGHI: " << sipHashRef("ABCDEFGHI", 9) << std::endl;
	std::cout << "ABCDEFGHI: " << SipHash24runtime("ABCDEFGHI") << std::endl;
	std::cout << "ABCDEFGHI: " << SipHash24("ABCDEFGHI") << std::endl;
	std::cout << "ABCDEFGHI: " << as_integer(Bla::ABCDEFGHI) << std::endl;

	std::cout << "ABCDEFGHIJ: " << sipHashRef("ABCDEFGHIJ", 10) << std::endl;
	std::cout << "ABCDEFGHIJ: " << SipHash24runtime("ABCDEFGHIJ" ) << std::endl;
	std::cout << "ABCDEFGHIJ: " << SipHash24("ABCDEFGHIJ") << std::endl;
	std::cout << "ABCDEFGHIJ: " << as_integer(Bla::ABCDEFGHIJ) << std::endl;

	std::cout << "ABCDEFGHIJK: " << sipHashRef("ABCDEFGHIJK", 11) << std::endl;
	std::cout << "ABCDEFGHIJK: " << SipHash24runtime("ABCDEFGHIJK" ) << std::endl;
	std::cout << "ABCDEFGHIJK: " << SipHash24("ABCDEFGHIJK") << std::endl;
	std::cout << "ABCDEFGHIJK: " << as_integer(Bla::ABCDEFGHIJK) << std::endl;

	std::cout << "ABCDEFGHIJKL: " << sipHashRef("ABCDEFGHIJKL", 12) << std::endl;
	std::cout << "ABCDEFGHIJKL: " << SipHash24runtime("ABCDEFGHIJKL" ) << std::endl;
	std::cout << "ABCDEFGHIJKL: " << SipHash24("ABCDEFGHIJKL") << std::endl;
	std::cout << "ABCDEFGHIJKL: " << as_integer(Bla::ABCDEFGHIJKL) << std::endl;

	std::cout << "ABCDEFGHIJKLM: " << sipHashRef("ABCDEFGHIJKLM", 13) << std::endl;
	std::cout << "ABCDEFGHIJKLM: " << SipHash24runtime("ABCDEFGHIJKLM" ) << std::endl;
	std::cout << "ABCDEFGHIJKLM: " << SipHash24("ABCDEFGHIJKLM") << std::endl;
	std::cout << "ABCDEFGHIJKLM: " << as_integer(Bla::ABCDEFGHIJKLM) << std::endl;

	std::cout << "ABCDEFGHIJKLMN: " << sipHashRef("ABCDEFGHIJKLMN", 14) << std::endl;
	std::cout << "ABCDEFGHIJKLMN: " << SipHash24runtime("ABCDEFGHIJKLMN") << std::endl;
	std::cout << "ABCDEFGHIJKLMN: " << SipHash24("ABCDEFGHIJKLMN") << std::endl;
	std::cout << "ABCDEFGHIJKLMN: " << as_integer(Bla::ABCDEFGHIJKLMN) << std::endl;

	std::cout << "ABCDEFGHIJKLMNO: " << sipHashRef("ABCDEFGHIJKLMNO", 15) << std::endl;
	std::cout << "ABCDEFGHIJKLMNO: " << SipHash24runtime("ABCDEFGHIJKLMNO") << std::endl;
	std::cout << "ABCDEFGHIJKLMNO: " << SipHash24("ABCDEFGHIJKLMNO") << std::endl;
	std::cout << "ABCDEFGHIJKLMNO: " << as_integer(Bla::ABCDEFGHIJKLMNO) << std::endl;

	std::cout << "ABCDEFGHIJKLMNOP: " << sipHashRef("ABCDEFGHIJKLMNOP", 16) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOP: " << SipHash24runtime("ABCDEFGHIJKLMNOP") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOP: " << SipHash24("ABCDEFGHIJKLMNOP") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOP: " << as_integer(Bla::ABCDEFGHIJKLMNOP) << std::endl;

	std::cout << "ABCDEFGHIJKLMNOPQ: " << sipHashRef("ABCDEFGHIJKLMNOPQ", 17) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQ: " << SipHash24runtime("ABCDEFGHIJKLMNOPQ" ) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQ: " << SipHash24("ABCDEFGHIJKLMNOPQ") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQ: " << as_integer(Bla::ABCDEFGHIJKLMNOPQ) << std::endl;

	std::cout << "ABCDEFGHIJKLMNOPQR: " << sipHashRef("ABCDEFGHIJKLMNOPQR", 18) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQR: " << SipHash24runtime("ABCDEFGHIJKLMNOPQR") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQR: " << SipHash24("ABCDEFGHIJKLMNOPQR") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQR: " << as_integer(Bla::ABCDEFGHIJKLMNOPQR) << std::endl;

	std::cout << "ABCDEFGHIJKLMNOPQRS: " << sipHashRef("ABCDEFGHIJKLMNOPQRS", 19) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRS: " << SipHash24runtime("ABCDEFGHIJKLMNOPQRS") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRS: " << SipHash24("ABCDEFGHIJKLMNOPQRS") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRS: " << as_integer(Bla::ABCDEFGHIJKLMNOPQRS) << std::endl;

	std::cout << "ABCDEFGHIJKLMNOPQRST: " << sipHashRef("ABCDEFGHIJKLMNOPQRST", 20) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRST: " << SipHash24runtime("ABCDEFGHIJKLMNOPQRST") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRST: " << SipHash24("ABCDEFGHIJKLMNOPQRST") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRST: " << as_integer(Bla::ABCDEFGHIJKLMNOPQRST) << std::endl;

	std::cout << "ABCDEFGHIJKLMNOPQRSTU: " << sipHashRef("ABCDEFGHIJKLMNOPQRSTU", 21) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTU: " << SipHash24runtime("ABCDEFGHIJKLMNOPQRSTU" ) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTU: " << SipHash24("ABCDEFGHIJKLMNOPQRSTU") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTU: " << as_integer(Bla::ABCDEFGHIJKLMNOPQRSTU) << std::endl;

	std::cout << "ABCDEFGHIJKLMNOPQRSTUV: " << sipHashRef("ABCDEFGHIJKLMNOPQRSTUV", 22) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUV: " << SipHash24runtime("ABCDEFGHIJKLMNOPQRSTUV" ) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUV: " << SipHash24("ABCDEFGHIJKLMNOPQRSTUV") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUV: " << as_integer(Bla::ABCDEFGHIJKLMNOPQRSTUV) << std::endl;

	std::cout << "ABCDEFGHIJKLMNOPQRSTUVW: " << sipHashRef("ABCDEFGHIJKLMNOPQRSTUVW", 23) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVW: " << SipHash24runtime("ABCDEFGHIJKLMNOPQRSTUVW") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVW: " << SipHash24("ABCDEFGHIJKLMNOPQRSTUVW") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVW: " << as_integer(Bla::ABCDEFGHIJKLMNOPQRSTUVW) << std::endl;

	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWX: " << sipHashRef("ABCDEFGHIJKLMNOPQRSTUVWX", 24) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWX: " << SipHash24runtime("ABCDEFGHIJKLMNOPQRSTUVWX" ) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWX: " << SipHash24("ABCDEFGHIJKLMNOPQRSTUVWX") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWX: " << as_integer(Bla::ABCDEFGHIJKLMNOPQRSTUVWX) << std::endl;

	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWXY: " << sipHashRef("ABCDEFGHIJKLMNOPQRSTUVWXY", 25) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWXY: " << SipHash24runtime("ABCDEFGHIJKLMNOPQRSTUVWXY") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWXY: " << SipHash24("ABCDEFGHIJKLMNOPQRSTUVWXY") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWXY: " << as_integer(Bla::ABCDEFGHIJKLMNOPQRSTUVWXY) << std::endl;

	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWXYZ: " << sipHashRef("ABCDEFGHIJKLMNOPQRSTUVWXYZ", 26) << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWXYZ: " << SipHash24runtime("ABCDEFGHIJKLMNOPQRSTUVWXYZ") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWXYZ: " << SipHash24("ABCDEFGHIJKLMNOPQRSTUVWXYZ") << std::endl;
	std::cout << "ABCDEFGHIJKLMNOPQRSTUVWXYZ: " << as_integer(Bla::ABCDEFGHIJKLMNOPQRSTUVWXYZ) << std::endl;

	return 0;
}
