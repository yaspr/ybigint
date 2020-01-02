#include <stdio.h>
#include <stdlib.h>

#include "bigint.h"

//
int main(int argc, char **argv)
{
  bigint_t *a = bigint_init_s("6700417");
  u8 *prime_msg[] = { "Composite", "Prime" };
  
  u8 p = bigint_is_probably_prime(a);

  printf("%s\n", prime_msg[p]);
  
  return 0;
  
  bigint_lib_init(0);
  
  //Generate 50 & 49 digits long numbers (random)
  //bigint_t *a = bigint_init_r(50);
  bigint_t *b = bigint_init_r(50);
  
  //Create an empty number
  bigint_t *c = bigint_init_z(0);
  
  //Add
  bigint_add_u(a, b, c);

  bigint_print(a);
  printf(" + ");
  bigint_print(b);
  printf(" = ");
  bigint_print(c);
  printf("\n\n");
  
  //Subtract
  bigint_sub_i(a, b, c);

  bigint_print(a);
  printf(" - ");
  bigint_print(b);
  printf(" = ");
  bigint_print(c);
  printf("\n\n");

  //Multiply
  bigint_mul_u(a, b, c);

  bigint_print(a);
  printf(" * ");
  bigint_print(b);
  printf(" = ");
  bigint_print(c);
  printf("\n\n");

  //Compare & divide

  u8 cmp;
  i64 pos;
  bigint_t *q = bigint_init_z(0);
  bigint_t *r = bigint_init_z(0);
  
  //If a > b
  if ((cmp = bigint_compare_i(a, b, &pos)) == 1)
    {
      printf("a / b\n");
      bigint_div_i(a, b, q, r);
      
      bigint_print(a);
      printf(" / ");
      bigint_print(b);
      printf(" => q: ");
      bigint_print(q);
      printf(" r: ");
      bigint_print(r);
      printf("\n\n");
    }
  else //If a < b
    if (cmp == 2)
      {
	printf("b / a\n");

	bigint_div_i(b, a, q, r);
      
	bigint_print(b);
	printf(" / ");
	bigint_print(a);
	printf(" => q: ");
	bigint_print(q);
	printf(" r: ");
	bigint_print(r);
	printf("\n\n");
      }
  
  //Factorial 
  bigint_t *d = bigint_init_s("500");

  bigint_fact_u(d, c);
  
  bigint_print(d);
  printf("! = ");
  bigint_print(c);
  printf("\n\n");

  bigint_free(&d);

  d = bigint_init_u(2);
  
  //Power
  u64 _e = 35;
  bigint_t *e = bigint_init_u(35);
  bigint_t *f = bigint_init_u(0);
  
  bigint_print(d);
  printf("^");
  bigint_print(e);

  bigint_pow_i(d, e, f);
  
  printf(" = ");
  bigint_print_meta(f);
  printf("\n\n");

  //Power
  bigint_print(d);
  printf("^%llu", _e);

  //
  bigint_pow_n(d, _e, f);
  
  printf(" = ");
  bigint_print_meta(f);
  printf("\n\n");

  //Freeing memory
  bigint_free(&a);
  bigint_free(&b);
  bigint_free(&c);
  bigint_free(&d);

  bigint_free(&q);
  bigint_free(&r);

  bigint_free(&e);

  //Fibonacci sequence
  a = bigint_init_u(0);
  b = bigint_init_u(1);
  c = bigint_init_z(0);

  //
  FILE *fp = fopen("fibo1000.txt", "wb");

  if (!fp)
    return printf("Cannot create Fibonacci file!\n"), 1;

  fprintf(fp, "0\t\t 0\n1\t\t 1\n");
  
  //
  for (u64 i = 2; i < 1000; i++)
    {
      printf("%llu\t\t ", i);
      
      bigint_add_u(a, b, c);

      (*a) = (*b);
      (*b) = (*c);
      
      bigint_print(c);
      putchar('\n');

      fprintf(fp, "%llu\t\t ", i);
      bigint_print_file(c, fp);
      fprintf(fp, "\n");
    }

  fclose(fp);
  
  return 0;
}
