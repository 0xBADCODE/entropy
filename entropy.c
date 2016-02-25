/* Entropy- Determine randomness / file analysis
 * testfileTRNG and testfilePRNG
 *
 * gcc entropy.c -o entropy -ggdb -Wall -lm
 *
 * TODO
 * FIX PI ESTIMATE
 * Copyright (c) 2016 Thomas Hand <th6045@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define byte sizeof(unsigned char)
#define BUFFER_SIZE 2

unsigned int i = 0, k = 0, m = 0, total_bytes = 0, total_bset = 0, alpha_count = 0;
static unsigned int freq[256], bset[256];
double bitpcent = 0, amv = 0, pi = 0, chi_sq = 0, lang_chi_sq= 0, entropy = 0, diff = 0;

static double p[256], p_alpha[26], p_english[26] = { 	
										8.167, 1.492, 2.782, 4.253, 12.702, 2.228, 2.015, 6.094, 6.966, 0.153, 0.772, 4.025, 2.406,
										6.749, 7.507, 1.929, 0.095, 5.987, 6.327, 9.056, 2.758, 0.978, 2.361, 0.150, 1.974, 0.074 };

unsigned char b[BUFFER_SIZE] = {'\0'}, alphabet[52] = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 
														'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
														'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 
														'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z' };
FILE *f = NULL;

/* log2 function: log2(x) = logy(x) / logy(2) */
double log2(double x)
{
	return (log10(x) / log10(2.0));
}

/* Monte Carlo estimate for pi */
void monte_carlo(double x, double y)
{
	x /= 256.0;
	y /= 256.0;

	/* Does (x, y) satisfy circle equation? */
	if((pow(x, 2.0) + pow(y, 2.0)) <= 1) m++;

	//printf("[ENTROPY]: x = %f, y = %f, m = %d, pi = %f\n", x, y, m, 4 * m / ((double)total_bytes/2)); //DEBUG
	//if(4.0 * m / ((double)total_bytes/2) < 3.14159265) //exit(-1); //DEBUG
	return;
}

double chi_squared(double obs, double expt)
{
	return pow((obs - expt), 2.0) / expt;
}

int check_array(int x, unsigned char *array){
    int j;
    for(j = 0; j < sizeof array; j++) {
        if (array[j] == x)
            return 1;
    }
    return 0;
}

int main(int argc, char *argv[])
{
	printf("[ENTROPY]: Analysing file...\n\n");

	if(argv[1] != NULL) 
		if (!(f = fopen(argv[1], "rb")))
			return -1;
		else ;
	else
		f = stdin;

	/* Read bytes into memory and process, repeat until no more bytes */
	while(fread(b, byte, BUFFER_SIZE, f) != 0) 
	{
		for(i = 0; i < BUFFER_SIZE; i++)
		{
			total_bytes++;

			/* Categorise distribution */
			freq[b[i]]++;

			/* Count only alphabet bytes */
			alpha_count += check_array(b[i], alphabet);
		}

		monte_carlo(b[0], b[1]);
		memset(&b, 0, BUFFER_SIZE);
	}

	for(i = 0; i < 256; i++)
	{
		/* Table of bits set per byte 0-255 */
		k = i;
		for(bset[i] = 0; k > 0; k >>= 1)
		{
			bset[i] += k & 0x1;
		}
		//printf("[ENTROPY]: Value[%d]: bits set = %d\n", i, bset[i]); //DEBUG

		/* Total bits set */
		total_bset += freq[i] * bset[i];

		/* Calculate arithmetic mean value */
		amv += i * freq[i];

		/* Chi squared distribution */
		chi_sq += chi_squared(freq[i], total_bytes / 256.0);

		//printf("[ENTROPY]: Chi-Squared Value: %f\n", chi_sq); //DEBUG

		/* Calculate entropy */
		p[i] = ((double)freq[i]) / total_bytes;
		if (p[i] > 0)
		{
			entropy += p[i] * log2(1 / p[i]);
		}
		//printf("[ENTROPY]: Probability: %f, Entropy: %f\n", p[i], entropy); //DEBUG
	}

	printf("Number of samples = %d\nNumber of bits set = %d", total_bytes, total_bset);
	
	/* Calculate 50/50 bits */
	bitpcent = total_bset / (total_bytes * 0.08);
	printf(" (%.4f/%.4f)%% (50/50 = random)\n", bitpcent, 100 - bitpcent);
	printf("Bit Ratio = %.8f (Ratio of 1.0 = random)\n", total_bset / (total_bytes * 8.0 - total_bset));

	printf("Arithmetic Mean Value = %0.3f (127.5 = random)\n", amv / total_bytes);

	pi = 4 * ((double)m / (total_bytes / 2)); // m/n = pi/4
	printf("Monte Carlo estimate of Pi = %0.8f (%.2f%% error) (TODO: FIX ERROR)\n", pi, 100 * (fabs((M_PI - pi) / M_PI)));

	printf("Chi-Squared distribution = %.3f (lower is more random)\n", chi_sq);

	printf("Entropy = %.6f bits per byte (8.0 = random) and could be further compressed by %.1f%%\n", entropy, (8 - entropy)/0.08);

	/* Plot byte frequency distribution */
	printf("\n\t--Frequency Analysis--\n");
/*	for (i = 0; i < 256; i++) //get largest freq element k
	{
		if (freq[i] > k)
		{
			k = freq[i];
		}
	}

	for(i = 0; i < 256; i++)
	{
		printf("0x%02x:", i);

		for(c = 0; c < (100 * freq[i] / k); c++)
		{
			printf("#");
		}

		printf("\n");
	}
*/
	/* Do frequency analysis */
    /* Calculate observed language probabilities and chi_sq */
	printf("\n\tSample:\t\tEnglish:\t%% diff:\n");

	for(i = 0; i < 26; i++)
	{
		p_english[i] *= 0.01; //normalise to 1
		p_alpha[i] = (((double)freq['a' + i]) + (double)freq['A' + i]) / alpha_count;
		lang_chi_sq += chi_squared((double)freq['a' + i], p_english[i] * alpha_count);
		diff = (fabs(p_alpha[i] - p_english[i]) / p_english[i]);

		printf("[%c]\t%.5f\t\t%.5f\t\t%.2f\n", 'a' + i, p_alpha[i], p_english[i], 100 * diff);
	}

	//printf("[ENTROPY]: CHI-SQUARED TEST ON LANGUAGE: %.3f\n", lang_chi_sq); //DEBUG

	/* Arbitrary tests */
	if(chi_sq > 500000 && lang_chi_sq < 20000)
		printf("\nThis file contains text and is not encrypted or random\n\n");
	else
		if(chi_sq > 30000)
			printf("\nThis file (probably) contains text and is not encrypted or random\n\n");
	else
		if(chi_sq > 1000)
			printf("\nThis file is not likely to be encrypted or random\n\n");
	else
		if(chi_sq > 300)
			printf("\nThis file is (probably) encrypted or sufficiently compressed to appear random\n\n");
	else
		printf("\nThis file is encrypted or random\n\n");

    return 0;
}
