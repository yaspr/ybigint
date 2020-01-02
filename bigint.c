//

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "bigint.h"

//
u8 randxy(u8 x, u8 y)
{ return x + (rand() % (y - x + 1)); }

//
void bigint_lib_init(u64 s)
{
  if (s && s != 1)
    {
      srand(s);
    }
  else
    {
      srand(getpid());
    }
}

//Init from: unsigned & signed integers, string, and file

//Create a bigint of n digits initialized to zero
bigint_t *bigint_init_z(u64 n)
{
  bigint_t *b = malloc(sizeof(bigint_t));

  if (b)
    {
      b->sign = 0;
      b->ndigits = n;

      if (n)
	memset(b->digits, 0, sizeof(u8) * n);
      
      return b;
    }
  else
    return NULL;
}

//Init from unsigned
bigint_t *bigint_init_u(u64 v)
{
  u64 i = 0;
  bigint_t *b = malloc(sizeof(bigint_t));

  if (b)
    {
      while (v)
	{
	  b->digits[i++] = v % 10;
	  v /= 10;
	}
      
      b->sign = 0;
      b->ndigits = i;
      
      return b;
    }
  else
    return NULL;
}

//Init from signed
bigint_t *bigint_init_i(i64 v)
{
  bigint_t *b = bigint_init_u(BIGINT_MACRO_ABS(v));
  
  if (b)
    b->sign = (v < 0);
  
  return b;
}

//Init from string with handling of first caharcter as a sign, if present
bigint_t *bigint_init_s(u8 *s)
{
  u8 *p = s;
  
  if (p)
    {
      bigint_t *b = malloc(sizeof(bigint_t));
      
      if (b)
	{
	  //Handle negative sign
	  if (p[0] == '-')
	    {
	      b->sign = 1;
	      p++;
	    }
	  else //Handle positive sign & non specified sign
	    {
	      b->sign = 0;
	      
	      //Move only if the character is +
	      p += (p[0] == '+');
	    }

	  //Compute the length after making sure sign characters have been handled
	  u64 len = strlen(p);

	  //
	  b->ndigits = len;
	  
	  //Copy in reverse an convert from digit character to digit values ('0' ==> 0, '1' ==> 1, ...) 
	  for (i64 i = len - 1, j = 0; i >= 0; i--, j++)
	    b->digits[j] = p[i] - '0';
	  
	  return b;
	}
      else
	return NULL;
    }
  else
    return NULL;
}

//Init from file
bigint_t *bigint_init_f(u8 *fname)
{
  FILE *fp = fopen(fname, "rb");

  if (fp)
    {
      u8 tmp[BIGINT_MAX_DIGITS + 1];
      u64 n = fread(tmp, sizeof(u8), BIGINT_MAX_DIGITS, fp);

      fclose(fp);
      
      bigint_t *b = bigint_init_s(tmp);

      return b;
    }
  else
    return NULL;
}

//Init with pseudo random values
bigint_t *bigint_init_r(u64 n)
{
  if (n)
    {
      bigint_t *b = malloc(sizeof(bigint_t));

      if (b)
	{
	  b->sign = 0;
	  b->ndigits = n;
	  
	  for (u64 i = 0; i < n - 1; i++)
	    b->digits[i] = randxy(0, 9);
	  
	  b->digits[n - 1] = randxy(1, 9);
	  
	  return b;
	}
      else
	return NULL;
    }
  else
    return NULL;
}

//Comparison

//Check if a bigint is 0
u8 bigint_is_zero(bigint_t *a)
{
  //If number of digits is 0 then assume empty
  if (a->ndigits)
    {
      u64 n = 0;
      
      while (n < a->ndigits && !a->digits[n])
	n++;
      
      return (n == a->ndigits);
    }
  else
    return 1;
}

//
u8 bigint_is_even(bigint_t *a)
{
  return !(a->digits[0] & 1);
}

//
u8 bigint_is_odd(bigint_t *a)
{
  return (a->digits[0] & 1);
}

//Check if two bigint values are equal (pos is the index of the mismatching digits) 
u8 bigint_equal_u(bigint_t *a, bigint_t *b, i64 *pos)
{
  (*pos) = -1;
  
  //If sizes aren't equal then values cannot be equal
  if ((a->ndigits > b->ndigits) || (a->ndigits < b->ndigits))
    return 0;
  else
    {
      //
      i64 i = 0;
      u8 same = 1;
      u64 n = a->ndigits; //or n = b->ndigits;
      
      //Using a xor to avoid using a compare
      //same = !(a xor b)
      //if a xor b == 0 then same = 1
      //If a xor b != 0 then same = 0
      while (same && i < n)
	{
	  same = !(a->digits[i] ^ b->digits[i]);
	  i += same;
	}
      
      (*pos) = i;
      
      return same;
    }
}

//
u8 bigint_equal_i(bigint_t *a, bigint_t *b, i64 *pos)
{
  *(pos) = -1;

  //If one is positive and the other negative ==> not equal
  if ((a->sign > b->sign) || (a->sign < b->sign)) 
    return 0;
  else
    return bigint_equal_u(a, b, pos);
}

//Check if two bigint values are 
u8 bigint_compare_u(bigint_t *a, bigint_t *b, i64 *pos)
{
  (*pos) = -1;
  
  //Check by size
  if (a->ndigits > b->ndigits)
    return 1; //First operand is larger than the second
  else
    if (a->ndigits < b->ndigits)
      return 2; //Second operand is larger than the first 
    else
      {
	//If size is similar

	//Check equality
	if (bigint_equal_u(a, b, pos))
	  return 0;
	else
	  {
	    u8 same = 1;
	    u64 n = a->ndigits;
	    i64 i = n - 1;

	    //
	    while (same)
	      {
		same = !(a->digits[i] ^ b->digits[i]);
		i -= same;
	      }
	    
	    (*pos) = i;
	    
	    return (a->digits[i] > b->digits[i]) ? 1 : 2;
	  }
      }
}

//
u8 bigint_compare_i(bigint_t *a, bigint_t *b, i64 *pos)
{
  (*pos) = -1;

  if (a->sign == 0 && b->sign == 1) //a positive & b negative
    return 1;
  else
    if (a->sign == 1 && b->sign == 0) //a negative & b positive
      return 2;
    else
      {
	u8 cmp = bigint_compare_u(a, b, pos);

	if (a->sign && b->sign)
	  {
	    if (cmp == 1)
	      return 2;
	    else
	      if (cmp == 2)
		return 1;
	      else
		return cmp;
	  }
	else
	  return cmp;
      }
}

//Addition

//Addition of two digits with carry
static inline u8 _add_(u8 a, u8 b, u8 *c)
{
  u8 r = a + b + (*c);

  (*c) = r / 10;

  return r % 10;
}

//Addition of two bigints with carry propagation: r = a + b
void bigint_add_u(bigint_t *a, bigint_t *b, bigint_t *r)
{
  u8 *p = NULL;
  u8 carry = 0;
  u64 min_ndigits;
  u64 max_ndigits;

  //
  r->sign = 0;
  
  //
  if (a->ndigits < b->ndigits)
    {
      min_ndigits = a->ndigits;
      max_ndigits = b->ndigits;
      p = b->digits;
    }
  else
    {
      min_ndigits = b->ndigits;
      max_ndigits = a->ndigits;
      p = a->digits;
    }
  
  //
  for (u64 i = 0; i < min_ndigits; i++)
    r->digits[i] = _add_(a->digits[i], b->digits[i], &carry);
  
  //
  for (u64 i = min_ndigits; i < max_ndigits; i++)
    r->digits[i] = _add_(p[i], 0, &carry);
  
  //
  if (carry)
    {
      r->digits[max_ndigits] = carry;
      r->ndigits = max_ndigits + 1;
    }
  else
    r->ndigits = max_ndigits;
}

//Signed addition
// (+, +): a + b =  |a| +  |b|; 
// (-, -): a + b = -|a| -  |b| = -(|a| + |b|);
// (-, +): a + b = -|a| +  |b| =  (|b| - |a|);
// (+, -): a + b = +|a| -  |b| =   |a| - |b|;
void bigint_add_i(bigint_t *a, bigint_t *b, bigint_t *r)
{
  if (bigint_is_zero(a))
    (*r) = (*b);
  else
    if (bigint_is_zero(b))
      (*r) = (*a);
    else
      {
	u8 s_a = a->sign;
	u8 s_b = b->sign;
	
	if (s_a == 0 && s_b == 0)
	  bigint_add_u(a, b, r);
	else
	  if (s_a == 1 && s_b == 1)
	    {
	      bigint_add_u(a, b, r);
	      r->sign = 1;
	    }
	  else
	    if (s_a == 1 && s_b == 0)
	      {
		a->sign = 0;
		bigint_sub_i(b, a, r);
		a->sign = 1;
	      }
	    else
	      if (s_a == 0 && s_b == 1)
		{
		  b->sign = 0;
		  bigint_sub_i(a, b, r);
		  b->sign = 1;
		}
      }
}

//
static inline u8 _sub_(u8 a, u8 b, u8 *c)
{
  i8 r = (a - (b + (*c)));
  u8 s = (r < 0);

  //If a < b then compensate with 10
  r = r + ((s) ? (10) : (0));

  //If compensation applied, set carry for the next iteration
  (*c) = s;

  return r;
}

//Subtraction - this function assumes (a >= b)
void bigint_sub_u(bigint_t *a, bigint_t *b, bigint_t *r)
{
  u8 *p = NULL;
  u8 carry = 0;
  u64 min_ndigits;
  u64 max_ndigits;

  //
  r->sign = 0;
  
  //
  if (a->ndigits < b->ndigits)
    {
      min_ndigits = a->ndigits;
      max_ndigits = b->ndigits;
      p = b->digits;
    }
  else
    {
      min_ndigits = b->ndigits;
      max_ndigits = a->ndigits;
      p = a->digits;
    }
  
  //
  for (u64 i = 0; i < min_ndigits; i++)
    r->digits[i] = _sub_(a->digits[i], b->digits[i], &carry);
  
  //
  for (u64 i = min_ndigits; i < max_ndigits; i++)
    r->digits[i] = _sub_(p[i], 0, &carry);
  
  //
  i64 n = (i64)max_ndigits - 1;
  
  //Remove trailing zeroes
  while (n >= 0 && !r->digits[n])
    n--;

  //When n == -1 ==> a - b = 0
  if (n < 0)
    r->ndigits = 1;
  else
    r->ndigits = (u64)n + 1;  
}

//if (|a| > |b|) ==> r = a - b
//if (|a| < |b|) ==> r = a - b = -(b - a) ==> calculate (b - a) and set r sign to '-'
void bigint_compare_swap_sub_u(bigint_t *a, bigint_t *b, bigint_t *r)
{
  i64 pos;
  u8 cmp = bigint_compare_u(a, b, &pos); //(compare)

  if (cmp == 1)   //a > b ==> a - b
    bigint_sub_u(a, b, r); //(sub)
  else
    if (cmp == 2) //a < b ==> a - b = -(b - a) (swap)
      {
	bigint_sub_u(b, a, r); //(sub)
	r->sign = 1; //
      }
    else
      if (cmp == 0)
	{
	  r->sign = 0;
	  r->ndigits = 1;
	  r->digits[0] = 0;
	}
}

//Signed subtraction
// (+, +): a - b = +|a| - +|b|; 
// (-, -): a - b = -|a| - -|b| = (-|a| + |b|) = |b| - |a|;
// (-, +): a - b = -|a| - +|b| = -(|a| + |b|);
// (+, -): a - b = +|a| - -|b| = |a| + |b|;
void bigint_sub_i(bigint_t *a, bigint_t *b, bigint_t *r)
{
  //0 - b ==> b->sign = !b->sign; (invert the sign)
  if (bigint_is_zero(a))
    {
      (*r) = (*b);
      r->sign = !r->sign;
    }
  else //a - 0
    if (bigint_is_zero(b))
      {
	(*r) = (*a);
      }
    else
      {
	u8 s_a = a->sign;
	u8 s_b = b->sign;
	
	if (s_a == 0 && s_b == 0)   //(+, +) 
	  bigint_compare_swap_sub_u(a, b, r);
	else
	  if (s_a == 1 && s_b == 1) //(-, -)
	    bigint_compare_swap_sub_u(b, a, r);
	  else
	    if (s_a == 1 && s_b == 0) //(-, +)
	      {
		bigint_add_u(a, b, r);
		r->sign = 1;
	      }
	    else
	      bigint_add_u(a, b, r);
      }
}

//Multiplication

//Unsigned multiply
void bigint_mul_2(bigint_t *a, bigint_t *r)
{
  if (bigint_is_zero(a))
    {
      r->sign = 0;
      r->ndigits = 1;
      r->digits[0] = 0;
    }
  else
    {
      u8 _r, _m, _c = 0;
      u64 n = a->ndigits;
      
      r->sign = a->sign;
      r->ndigits = a->ndigits;
      
      //
      for (u64 i = 0; i < n; i++)
	{
	  _m = (a->digits[i] << 1) + _c;
	  
	  _r = (_m % 10);
	  _c = (_m / 10);
	  
	  r->digits[i] = _r;
	}
      
      if (_c)
	{
	  r->digits[r->ndigits] = _c;
	  r->ndigits++;
	}
    }
}

/*
  Ancient Egyptian Algorithm:
  ---------------------------
  
  Suppose we want to perform the following operation: 135 * 56

  First, convert 135 into its binary form:

  135 = 128 + 7 = 128 + 4 + 2 + 1 ==> 135 = 10000111
  
  Then, using the binary representation we perform the following:
  
  Bits | State | Weigths  | Final multiply   
  0    | 1     | 56 * 2^0 | 56 * 2^0
  1    | 1     | 56 * 2^1 | 56 * 2^0 + 56 * 2^1
  2    | 1     | 56 * 2^2 | 56 * 2^0 + 56 * 2^1 + 56 * 2^2
  3    | 0     | skip     | 0
  4    | 0     | skip     | 0
  5    | 0     | skip     | 0
  6    | 0     | skip     | 0
  7    | 1     | 56 * 2^7 | 56 * 2^0 + 56 * 2^1 + 56 * 2^2 + 56 * 2^7
  
  56 * 2^0 + 56 * 2^1 + 56 * 2^2 + 56 * 2^7 = 56 + 112 + 224 + 7168 = 7560

  Mathematical notation:
  

          nbits - 1
           ___
          \
  a * b =  |    bin(a)[i] * b * 2^i
          /___
	  
          i = 0

  For the previous example, we can set: a = 135 & b = 56
  
  bin(a = 135) = 10000111
  The number of bits (nbits) = 8

  Bit position   7 6 5 4 3 2 1 0
  bin(a = 135)   1 0 0 0 0 1 1 1
  
  a * b = 1 * 56 * 2^0 + 1 * 56 * 2^1 + 1 * 56 * 2^2 + 0 * 56 * 2^3 + 
          0 * 56 * 2^4 + 0 * 56 * 2^5 + 0 * 56 * 2^6 + 1 * 56 * 2^7
	
	= 56 + 56 * 2 + 56 * 4 + 56 * 128 
	
	= 56 + 112 + 224 + 7168

	= 7560


  Karatsuba's Algorithm:
  -----------------------

  X denotes the base
  m denotes the base exponent
  
  a = a1 * X^m + a0
  b = b1 * X^m + b0
  
  i.e. 135 = 13 * 10^1 + 5 = 1 * 10^2 + 35
  
  a * b = (a1 * X^m + a0) * (b1 * X^m + b0)
  
        = (a1 * b1) * X^(2m) + (a1 * b0) * X^m + (a0 * b1) * X^m + (a0 * b0)
        = (a1 * b1) * X^(2m) + (a1 * b0 + a0 * b1) * X^m + (a0 * b0)
	= (a1 * b1) * X^(2m) + (a1 * b0 + a0 * b1) * X^m + (a0 * b0)
	
	with (a1 * b0 + a0 * b1) = (a0 + a1) * (b0 + b1) - (a1 * b1 + a0 * b0)
	
	(this allows to move from a [* + *] pattern to a [+ * + - +] with only one multiply - see [0])
	
  a * b = (a1 * b1) * X^(2m) + ((a0 + a1) * (b0 + b1) - (a1 * b1 + a0 * b0)) * X^m + (a0 * b0)
  
  Let:
           A = (a1 * b1), 
	   B = (a0 * b0), 
	   C = (a0 + a1), 
	   D = (b0 + b1), 
	   E = C * D, 
	   F = A + B, 
           G = E - F
	   
  a * b = A * X^(2m) + G * X^m + B

  
  Let K be a recursive Karatsuba multiplication function:
  
  a * b = K(a, b)
  
  K(a, b) = A * X^(2m) + G * X^m + B

  Then:

          A = K(a1, b1), 
	  B = K(a0, b0), 
	  C = (a0 + a1), 
	  D = (b0 + b1) 
	  E = K(C, D), 
	  F = A + B 
	  G = E - F

  Three calls to K are performed per one call to K. 

  [0] NOTE:
  ---------
  
  Suppose we stop at this step:
  
  a * b = (a1 * b1) * X^(2m) + (a1 * b0 + a0 * b1) * X^m + (a0 * b0)

  Let:
         A = (a1 * b1),
	 B = (a0 * b0),
	 C = (a1 * b0),
	 D = (a0 * b1),
	 E = C + D

  After applying the Karatsuba function K:

         A = K(a1, b1),
	 B = K(a0, b0),
	 C = K(a1, b0),
	 D = K(a0, b1),
	 E = C + D

  Here, four calls to K are performed per one call to K.  

  Grade-school mutiplication algorithm:
  -------------------------------------

  Let's multiply abcd by xyz where a, b, c, d, x, y, and z denote the digits.
  
  len(abcd) denotes the number of digits in a number. In this case len(abcd) = 4, 
  and len(xyz) = 3. 

                           abcd
   *                        xyz
    ---------------------------
              za + zb + zc + zd  level 0 \
         ya + yb + yc + yd +  0  level 1  | Perform horizontal add for each level 
    xa + xb + xc + xd +  0 +  0  level 2 /
    ---------------------------
    Vertical add for each column

    This algorithm requires len(abcd) * len(xyz) multiplications (4 * 3 = 12 in this example)
    and 12 horizontal additions.
    
    Another way to represent this multiplication is by using the following matrix:

    Multiply each row element by the columns elements and store the value at the (row, col)
    position.
    
      |   a |  b |  c |  d
    ----------------------
    x |  xa | xb | xc | xd --> Add horizontally and multiply by 100 (similar to previous level 2)
    y |  ya | yb | yc | yd --> Add horizontally and multiply by 10  (similar to previous level 1)
    z |  za | zb | zc | zd --> Add horizontally                     (similar to previous level 0)

    A better way is to perform the multiplication upside down:

      |   d |  c |  b |  a
    ----------------------
    z |  zd | zc | zb | za --> Add horizontally                     (level 0)
    y |  yd | yc | yb | ya --> Add horizontally and multiply by 10  (level 1)
    x |  xd | xc | xb | za --> Add horizontally and multiply by 100 (level 2)

    The result is (level 0) + (level 1) + (level 2). 
*/

//Unsigned multiply the ancient Egyptian way
void bigint_mul_u(bigint_t *a, bigint_t *b, bigint_t *r)
{
  u8 _r;
  bigint_t _a = (*a), _b = (*b);

  r->sign = 0;
  r->ndigits = 0;
  r->digits[0] = 0;
  
  //
  do
    {
      bigint_div_2(&_a, &_a, &_r);

      //Only accumulate when the bit is up
      if (_r)
	bigint_add_u(r, &_b, r);

      //Double the value (b = 2 * b or b = b + b)
      bigint_mul_2(&_b, &_b);
    }
  while (!bigint_is_zero(&_a));
}

//
void bigint_mul_i(bigint_t *a, bigint_t *b, bigint_t *r)
{
  u8 s_a = a->sign;
  u8 s_b = b->sign;
  
  bigint_mul_u(a, b, r);

  if (s_a == s_b)
    r->sign = 0;
  else
    r->sign = 1;
}

//Division

//Divide by 2
void bigint_div_2(bigint_t *a, bigint_t *q, u8 *r)
{
  if (bigint_is_zero(a))
    {
      q->sign = 0;
      q->ndigits = 1;
      q->digits[0] = 0;

      (*r) = 0;
    }
  else
    {
      i64 n = a->ndigits - 1;
      u8 _v = a->digits[n], _q, _r;
      
      q->sign = a->sign;
      q->ndigits = n + 1;
      
      while (n >= 0)
	{
	  _q = _v >> 1;
	  _r = _v & 1;
	  
	  q->digits[n] = _q;
	  
	  _v = (_r * 10) + a->digits[n - 1];
	  
	  n--;
	}
      
      (*r) = _r;
      
      //Remove trailing zeros
      if (a->ndigits > 1)
	{
	  n = a->ndigits - 1;
	  
	  while (n >= 0 && !q->digits[n])
	    n--;
	  
	  q->ndigits = n + 1;
	}
    }
}

//Brute force search for q - assumes a > b (VERY SLOW when a is VERY large and b VERY small)
void bigint_div_bf_i(bigint_t *a, bigint_t *b, bigint_t *q, bigint_t *r)
{
  if (bigint_is_zero(b))
    {
      printf("Bigint error: division by 0!\n");
      exit(1);
    }
  
  if (bigint_is_zero(a))
    {
      q->sign = 0;
      q->ndigits = 1;
      q->digits[0] = 0;

      r->sign = 0;
      r->ndigits = 1;
      r->digits[0] = 0;
    }
  else
    {
      //a / b ==> a = b * q + r
      u8 cmp;
      i64 pos;
      bigint_t *_m = bigint_init_u(0);
      bigint_t *_q = bigint_init_u(1);
      bigint_t *_1 = bigint_init_u(1); //Increment an inc and dec

      //m = q * b
      bigint_mul_i(b, _q, _m);
      
      //while (m < a)
      while ((cmp = bigint_compare_u(_m, a, &pos)) == 2)
      	{
	  //q++ ===> Very slow yet accurate convergence to the largest q
      	  bigint_add_i(_q, _1, _q);

	  //m = b * q
      	  bigint_mul_u(b, _q, _m);
      	}

      //if (m > a)
      if (cmp)
	{
	  //q-- ==> adjust to the proper value 
	  bigint_sub_i(_q, _1, _q);

	  //m = b * q
	  bigint_mul_i(b, _q, _m);

	  //r = a - m
	  bigint_sub_i(a, _m, r);
	}
      else //if (m == a) ==> a is divisible by b
	{
	  //Rest is 0
	  r->sign = 0;
	  r->ndigits = 1;
	  r->digits[0] = 0;
	}

      //Output q
      (*q) = (*_q);

      //+ ==> 0 and - ==> 1
      //if similar signs are multiplied   ==> +
      //if different signs are multiplied ==> -
      q->sign = (a->sign ^ b->sign);
      
      free(_m);
      free(_q);
      free(_1);
    }
}

//
void bigint_div_i(bigint_t *a, bigint_t *b, bigint_t *q, bigint_t *r)
{
  if (bigint_is_zero(b))
    {
      printf("Bigint error: division by 0!\n");
      exit(1);
    }
  
  if (bigint_is_zero(a))
    {
      q->sign = 0;
      q->ndigits = 1;
      q->digits[0] = 0;

      r->sign = 0;
      r->ndigits = 1;
      r->digits[0] = 0;
    }
  else
    {
      /*
	q  = 1; //Current q
	pq = 1; //Previous q

	//
	m = b * q;
	
	//
	while (m < a)
	{
  	  pq = q;
	  q *= 2;    //Doubling to converge 
	  m = b * q;
	}
	
	q = pq;
	m = b * q;
	
	while (m < a)
	{
	  q++;
	  m = b * q;
	}

	r = a - m;
	
	return q, r ;
      */
      
      //
      u8 cmp;
      i64 pos;
      bigint_t *_q  = bigint_init_u(1);
      bigint_t *_m  = bigint_init_u(0);
      bigint_t *_1  = bigint_init_u(1);
      bigint_t *_pq = bigint_init_u(0);

      bigint_mul_u(b, _q, _m);
      
      //
      while ((cmp = bigint_compare_u(_m, a, &pos)) == 2)
	{
	  (*_pq) = (*_q);
	  
	  bigint_mul_2(_q, _q);
	  bigint_mul_u(b, _q, _m);
	}

      //If a = b * q
      if (cmp == 0)
	{
	  r->sign = 0;
	  r->ndigits = 1;
	  r->digits[0] = 0;
	}
      else
	{
	  //
	  (*_q) = (*_pq);

	  bigint_mul_u(b, _q, _m);
	  
	  //
	  while ((cmp = bigint_compare_u(_m, a, &pos)) == 2)
	    {
	      (*_pq) = (*_q);
	      
	      bigint_add_u(_q, _1, _q);
	      
	      bigint_mul_u(b, _q, _m);
	    }

	  //if (m > a)
	  if (cmp)
	    {
	      //q-- ==> adjust to the proper value 
	      bigint_sub_i(_q, _1, _q);
	      
	      //m = b * q
	      bigint_mul_i(b, _q, _m);
	      
	      //r = a - m
	      bigint_sub_i(a, _m, r);
	    }
	  else //if (m == a) ==> a is divisible by b
	    {
	      //Rest is 0
	      r->sign = 0;
	      r->ndigits = 1;
	      r->digits[0] = 0;
	    }
	  
	  //Output q
	  (*q) = (*_q);
	  
	  //+ ==> 0 and - ==> 1
	  //if similar signs are multiplied   ==> +
	  //if different signs are multiplied ==> -
	  q->sign = (a->sign ^ b->sign);	  
	}
      
      //
      free(_q);
      free(_1);
      free(_m);
      free(_pq);
    }
}

//Power (a^p ==> multiply a by itself p - 1 times
void bigint_pow_n(bigint_t *a, u64 p, bigint_t *r)
{
  if (p)
    {
      (*r) = (*a);

      p--;

      while (p)
	{
	  bigint_mul_i(r, a, r);
	  p--;
	}
    }
  else
    {
      r->sign = 0;
      r->ndigits = 1;
      r->digits[0] = 1;
    }
}

//
void bigint_pow_i(bigint_t *a, bigint_t *p, bigint_t *r)
{
  if (!bigint_is_zero(p))
    {
      bigint_t *_1 = bigint_init_u(1);
      
      (*r) = (*a);

      bigint_sub_i(p, _1, p);

      while (!bigint_is_zero(p))
	{
	  bigint_mul_i(r, a, r);
	  bigint_sub_i(p, _1, p);
	}

      free(_1);
    }
  else
    {
      r->sign = 0;
      r->ndigits = 1;
      r->digits[0] = 1;
    }
}

/*
  Primality test - Miller-Rabin:
  ------------------------------
  
  n - 1 = m . 2^k;   (if k = 1, then n is composite)
    
  T[0] = 2^m mod n;  (if T0 = 1 or T = -1, then n is prime)
    
  for (i = 1; i < k -  1; k++)
  {
     T[i] = T[i - 1]^2 mod n;
     T[i] -= n;

     if (T[i] == -1)
        return prime;
     else
        if (T[i] == 1)
	   return composite
  }

  Very slow sometimes!
*/
u8 bigint_is_probably_prime(bigint_t *n)
{
  i8 _r;
  i64 pos;
  u64 k = 0;
  bigint_t *_1 = bigint_init_u(1);
  bigint_t *_2 = bigint_init_u(2);
  bigint_t *m = bigint_init_z(0);
  bigint_t *t = bigint_init_z(0);
  bigint_t *r = bigint_init_z(0);

  //
  bigint_sub_u(n, _1, m);

  //Calculate k & m
  while (bigint_is_even(m))
    {
      bigint_div_2(m, m, &_r);

      k++;
    }

  printf("m: ");
  bigint_print(m);
  putchar('\n');
  
  //2^m
  bigint_pow_i(_2, m, t);
  
  //t %= n  
  bigint_div_i(t, n, m, r);

  (*t) = (*r);

  //Not a prime
  if (k == 1)
    {
      printf("k: %llu\n", k);

      bigint_print(t);
      putchar('\n');
      
      free(_1);
      free(_2);
      free(m);
      free(t);
      free(r);
      
      return 0;;
    }

  for (u64 i = 1; !_r && i <= k - 1; i++)
    {
      bigint_mul_u(t, t, m);
      bigint_div_i(m, n, t, r);
      bigint_sub_i(r, n, t);

      //If t == 1 ==> n is composite
      if (!bigint_compare_i(t, _1, &pos))
	{
	  
	  free(_1);
	  free(_2);
	  free(m);
	  free(t);
	  free(r);
	  
	  return 0;
	}
      else
	{
	  //If t == -1 ==> n is prime
	  
	  _1->sign = 1;
	  
	  if (!bigint_compare_i(t, _1, &pos))
	    {
	      free(_1);
	      free(_2);
	      free(m);
	      free(t);
	      free(r);
	      
	      return 1;
	    }
	  
	  _1->sign = 0;
	}
    }
  
  return 0;
}

//Factorial - F(n) = n * (n - 1) * (n - 2) * 2
/*
  Can be parallelized.

  With 2 Threads
    
     Thread0: partial factorial f0 = n * (n - 1) * ... * ((n / 2) -  1)
     Thread1: partial factorial f1 = (n / 2) * ((n / 2) - 1) * ... * 2     
*/
void bigint_fact_u(bigint_t *a, bigint_t *r)
{
  i64 pos;
  bigint_t *_1 = bigint_init_u(1);
  bigint_t *_2 = bigint_init_u(2);
  bigint_t *r_1 = bigint_init_z(0);

  (*r) = (*a);
  (*r_1) = (*r);
  
  while (bigint_compare_u(r_1, _2, &pos))
    {
      bigint_sub_i(r_1, _1, r_1);
      bigint_mul_u(r, r_1, r);
    }
}

//
void bigint_free(bigint_t **b)
{
  free(*b);
  (*b) = NULL;
}

//Absolute value: |a| = +a
void bigint_abs_i(bigint_t *b)
{
  //Flip sign value (0 ==> 1, 1 ==> 0)
  b->sign = 0;
}

//Print sign only when negative 
void bigint_print(bigint_t *b)
{
  if (b->sign)
    printf("-");
  
  for (i64 i = b->ndigits - 1; i >= 0; i--)
    printf("%c", '0' + b->digits[i]);
}

//Print signed
void bigint_print_signed(bigint_t *b)
{
  if (b->sign)
    printf("-");
  else
    printf("+");
  
  for (i64 i = b->ndigits - 1; i >= 0; i--)
    printf("%c", '0' + b->digits[i]);
}

//Print meta data
void bigint_print_meta(bigint_t *b)
{
  printf("Number of digits    :\t %llu\n",  b->ndigits);
  printf("Size in Bytes       :\t %llu\n",  b->ndigits);
  printf("Size in KiBytes     :\t %llu\n", (b->ndigits >> 10));
  printf("Sign                :\t %c\n"  , (b->sign) ? '-' : '+');
  
  printf("Digits              :\t ");
  
  for (i64 i = b->ndigits - 1; i >= 0; i--)
    printf("%c", '0' + b->digits[i]);

  printf("\n");
}

//Write to file
u8 bigint_print_file(bigint_t *b, FILE *fp)
{  
  if (fp)
    {
      if (b->sign)
	fprintf(fp, "-");
      
      for (i64 i = b->ndigits - 1; i >= 0; i--)
	fprintf(fp, "%c", '0' + b->digits[i]);
      
      return 0;
    }
  else
    return 1;
}

