//
#define BIGINT_MAX_DIGITS 1024 * 500 //1/2 MiB

//
#define BIGINT_MACRO_ABS(v)     ((v) > 0)   ? (v) : (-(v))

//
#define BIGINT_MACRO_MIN(a, b)  ((a) < (b)) ? (a) : (b)
#define BIGINT_MACRO_MAX(a, b)  ((a) > (b)) ? (a) : (b)

//
typedef unsigned char       u8;
typedef unsigned short     u16;
typedef unsigned int       u32;
typedef unsigned long long u64;

//
typedef signed char       i8;
typedef signed short     i16;
typedef signed int       i32;
typedef signed long long i64;

//
typedef struct bigint_s {

  //Number of digits
  u64 ndigits;

  //Sign (0 ==> +, 1 ==> -)
  u8 sign;

  //String of digits in reverse order
  u8 digits[BIGINT_MAX_DIGITS];

} bigint_t;

//
void bigint_lib_init(u64 s);

//
bigint_t *bigint_init_z(u64 n);
bigint_t *bigint_init_u(u64 v);
bigint_t *bigint_init_i(i64 v);
bigint_t *bigint_init_r(u64 n);
bigint_t *bigint_init_s(u8 *s);
bigint_t *bigint_init_f(u8 *fname);

//
u8 bigint_is_zero(bigint_t *a);
u8 bigint_equal_u(bigint_t *a, bigint_t *b, i64 *pos);
u8 bigint_equal_i(bigint_t *a, bigint_t *b, i64 *pos);
u8 bigint_compare_u(bigint_t *a, bigint_t *b, i64 *pos);
u8 bigint_compare_i(bigint_t *a, bigint_t *b, i64 *pos);

//
void bigint_add_u(bigint_t *a, bigint_t *b, bigint_t *r);
void bigint_add_i(bigint_t *a, bigint_t *b, bigint_t *r);
void bigint_sub_u(bigint_t *a, bigint_t *b, bigint_t *r);
void bigint_sub_i(bigint_t *a, bigint_t *b, bigint_t *r);

//
void bigint_abs_i(bigint_t *b);

//
void bigint_mul_2(bigint_t *a, bigint_t *r);
void bigint_mul_u(bigint_t *a, bigint_t *b, bigint_t *r);
void bigint_mul_i(bigint_t *a, bigint_t *b, bigint_t *r);

//
void bigint_div_2(bigint_t *a, bigint_t *q, u8 *r);
void bigint_div_i(bigint_t *a, bigint_t *b, bigint_t *q, bigint_t *r);

//
void bigint_pow_n(bigint_t *a, u64 p, bigint_t *r);
void bigint_pow_i(bigint_t *a, bigint_t *p, bigint_t *r);

//
u8 bigint_is_probably_prime(bigint_t *n);

//
void bigint_fact_u(bigint_t *a, bigint_t *r);

//
void bigint_free(bigint_t **b);

//
void bigint_print(bigint_t *b);
void bigint_print_signed(bigint_t *b);
void bigint_print_meta(bigint_t *b);
u8 bigint_print_file(bigint_t *b, FILE *fp);



