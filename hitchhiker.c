#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <isa-l/erasure_code.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "hitchhiker.h"
#include "helper.h"
#include "base.h"

/*
文件填充后大小：2K*len/2 = K*len, len最小2，即最小 128B
对于2MB的文件，推荐len为：2MB/64=32
对于1GB的文件，推荐len为：1GB/64=16MB (16,777,216)
*/

u8* encode_matrix; // M*K, shared

int split_file(FILE *stream, int len, u8 **stripe_a, u8 **stripe_b)
{
    fseek(stream, 0, SEEK_END);
    int filesz = ftell(stream); // int 最大支持文件大小：2047MB
    if (K * len < filesz)
    {
        fprintf(stderr, "len is too small\n");
        return -1;
    }
    rewind(stream);

    u8 buf[K * len];
    memset(buf, 0, sizeof(buf));
    if (fread(buf, filesz, 1, stream) != 1)
    {
        fprintf(stderr, "fread error\n");
        return -1;
    }
    u8 *p = buf;
    for (int i = 0; i < K; i++)
    {
        memcpy(stripe_a[i], p, len / 2);
        p += len / 2;
        memcpy(stripe_b[i], p, len / 2);
        p += len / 2;
    }
    return 0;
}

/*
f1(a)	                    f1(b)
f2(0,0,0,a4,...,a10)^f2(b)	f2(b)^f2(a1,a2,a3,0,...,0)
f3(a)						f3(b)^f2(0,0,0,a4,a5,a6,0,...,0)
f4(a)						f4(b)^f2(0,...,0,a7,a8,a9,0)
*/

void get_part_stripe(u8 **stripe_a, u8 **part_a, const int *filter, int len)
{
    for (int i = 0; i < K; i++)
    {
        if (filter[i])
            memcpy(part_a[i], stripe_a[i], len);
        else
            memset(part_a[i], 0, len);
    }
}

void hhk_encode_parity(int len, int p, u8** stripe_a, u8** stripe_b, u8** out_a, u8** out_b) {
    int i;
    u8* part_a[K];
    u8* coded_a[1];
    for (i=0; i<K; i++) {
        part_a[i] = malloc(len);
    }
    coded_a[0]=malloc(len);

    if (p==K) {
        hhk_encode_1(len, K, stripe_a, out_a);
        hhk_encode_1(len, K, stripe_b, out_b);
    } else if (p==K+1) {
        hhk_encode_1(len, K+1, stripe_b, out_a);
        hhk_encode_1(len, K+1, stripe_b, out_b);
        get_part_stripe(stripe_a, part_a, hhk_f2_filter4, len);
        hhk_encode_1(len, K+1, part_a, coded_a);
        for (i=0; i<len; i++) {
            out_a[0][i] ^= coded_a[0][i];
        }
        get_part_stripe(stripe_a, part_a, hhk_f2_filter1, len);
        hhk_encode_1(len, K+1, part_a, coded_a);
        for (i=0; i<len; i++) {
            out_b[0][i] ^= coded_a[0][i];
        }
    } else if (p==K+2) {
        hhk_encode_1(len, K+2, stripe_a, out_a);
        hhk_encode_1(len, K+2, stripe_b, out_b);
        get_part_stripe(stripe_a, part_a, hhk_f2_filter2, len);
        hhk_encode_1(len, K+1, part_a, coded_a);
        for (i=0; i<len; i++) {
            out_b[0][i] ^= coded_a[0][i];
        }
    } else {
        hhk_encode_1(len, K+3, stripe_a, out_a);
        hhk_encode_1(len, K+3, stripe_b, out_b);
        get_part_stripe(stripe_a, part_a, hhk_f2_filter3, len);
        hhk_encode_1(len, K+1, part_a, coded_a);
        for (i=0; i<len; i++) {
            out_b[0][i] ^= coded_a[0][i];
        }
    }
}

void hhk_encode_data(int len, u8 **stripe_a, u8 **stripe_b, u8 **parity_a, u8 **parity_b)
{
    for (int i=0; i<P; i++) {
        hhk_encode_parity(len, i+K, stripe_a, stripe_b, &parity_a[i], &parity_b[i]);
    }
}

int hhk_decode_e1(int len, int e, int p, u8** substripe_left, u8* parity, u8** out) {
    int i, j;

    u8* sub_encode_matrix = malloc(K*K);
    u8* invert_matrix = malloc(K*K);
    u8* decode_matrix = malloc(1*K); // only select coefficients of index e, from invert_matrix
    u8* gtbls = malloc(1*K*32);
    j=0;
    for (i=0; i<K; i++) {
        if (i!=e) {
            memcpy(&sub_encode_matrix[K*(j++)], &encode_matrix[K*i], K);
        }
    }
    memcpy(&sub_encode_matrix[K*j], &encode_matrix[p*K], K);
    if (gf_invert_matrix(sub_encode_matrix, invert_matrix, K) != 0) {
        fprintf(stderr, "invert matrix error");
        return -1;
    }
    memcpy(decode_matrix, &invert_matrix[e*K], K);
    ec_init_tables(K, 1, decode_matrix, gtbls);

    u8* srcs[K];
    for (i=0; i<K-1; i++) srcs[i] = substripe_left[i];
    srcs[K-1] = parity;

    ec_encode_data(len, K, 1, gtbls, srcs, out);

    return 0;
}

void hhk_encode_1(int len, int p, u8** substripe, u8** out) {
    u8 gftbls[1*K*32];
    ec_init_tables(K, 1, &encode_matrix[p*K], gftbls);
    ec_encode_data(len, K, 1, gftbls, substripe, out);
}

void hhk_decode_a(int len, int e, u8** stripe_b, u8** stripe_a_left, u8* parity, u8** out) {
    int i, j, o;
    int p = K+1+e/3; // 0,1,2 K+1; 3,4,5 K+2; 6,7,8 K+3
    u8* fb[1];
    fb[0]=malloc(len);
    hhk_encode_1(len, p, stripe_b, fb);
    u8* fa = malloc(len);
    memcpy(fa, fb[0], len);
    for (i=0; i<len; i++) {
        fa[i] ^= parity[i];
    }

    const int* filter;
    if (e>=0 && e<=2) filter=hhk_f2_filter1;
    else if (e>=3 && e<=5) filter=hhk_f2_filter2;
    else filter=hhk_f2_filter3;

    u8* stripe_a_src[K-1];
    j=0;
    o=0;
    for (i=0; i<K; i++) {
        if (i==e) continue;
        else if (filter[i]) {
            stripe_a_src[j++]=stripe_a_left[o++];
        }
        else {
            stripe_a_src[j] = malloc(len);
            memset(stripe_a_src[j], 0, len);
            j++;
        }
    }

    hhk_decode_e1(len, e, K+1, stripe_a_src, fa, out); // fix f2 (K+1 in encode_matrix)
}

void hhk_decode_a_last(int len, u8** stripe_b, u8* parity2, u8* parity3, u8* parity4, u8** out) {
    int i;
    u8* f2b[1];
    u8* f3b[1];
    u8* f4b[1];
    f2b[0]=malloc(len);
    f3b[0]=malloc(len);
    f4b[0]=malloc(len);
    hhk_encode_1(len, K+1, stripe_b, f2b);
    hhk_encode_1(len, K+2, stripe_b, f3b);
    hhk_encode_1(len, K+3, stripe_b, f4b);
    for (i=0; i<len; i++) {
        f2b[0][i] ^= parity4[i];
        f3b[0][i] ^= parity2[i];
        f4b[0][i] ^= parity3[i];
    }

    for (i=0; i<len; i++) {
        f2b[0][i] ^= f3b[0][i];
        f2b[0][i] ^= f4b[0][i];
    }
    // now f2b[0]: f2(0,0,...0,a10)
    u8* srcs[K-1];
    for (i=0; i<K-1; i++) {
        srcs[i]=malloc(len);
        memset(srcs[i], 0, len);
    }
    hhk_decode_e1(len, K-1, K+1, srcs, f2b[0], out);
}

// simulate
void encode_file(char* filename) {
    struct stat st = {0};
    if (stat("a", &st)) {
        mkdir("a", 0775);
    }
    if (stat("b", &st)) {
        mkdir("b", 0775);
    }
    if (stat("pa", &st)) {
        mkdir("pa", 0775);
    }
    if (stat("pb", &st)) {
        mkdir("pb", 0775);
    }

    int i;

    FILE *stream = fopen(filename, "r");
    int full_len = 8; // 8*10=80B, 子块是 4B
    int len = full_len/2;

    u8 *stripe_a[K];
    u8 *stripe_b[K];
    for (i = 0; i < K; i++)
    {
        stripe_a[i] = malloc(len);
        stripe_b[i] = malloc(len);
    }
    split_file(stream, full_len, stripe_a, stripe_b);
    for (i = 0; i < K; i++)
    {
        printf("%d: ", i);
        printf("(%s) (%s)\n", stripe_a[i], stripe_b[i]);
    }

    u8 *parity_a[P];
    u8 *parity_b[P];
    for (i = 0; i < P; i++)
    {
        parity_a[i] = malloc(len);
        parity_b[i] = malloc(len);
    }
    hhk_encode_data(len, stripe_a, stripe_b, parity_a, parity_b);
    for (i=0; i<P; i++) {
        printf("%d: ", i);
        printf("(%s) (%s)\n", parity_a[i], parity_b[i]);
    }

    // write to file
    char name[10];
    FILE* os;
    for (i=0; i<K; i++) {
        sprintf(name, "a/a%d", i); // a0-a9
        os = fopen(name, "w");
        fwrite(stripe_a[i], len, 1, os);
        fclose(os); // BUG: if not close, fread error later

        sprintf(name, "b/b%d", i); 
        os = fopen(name, "w");
        fwrite(stripe_b[i], len, 1, os);
        fclose(os);
    }
    for (i=0; i<P; i++) {
        sprintf(name, "pa/pa%d", i); // pa0-pa3
        os = fopen(name, "w");
        fwrite(parity_a[i], len, 1, os);
        fclose(os);

        sprintf(name, "pb/pb%d", i);
        os = fopen(name, "w");
        fwrite(parity_b[i], len, 1, os);
        fclose(os);
    }
}

void get_block(int len, char* type, int index, u8* out) {
    char name[10];
    sprintf(name, "%s/%s%d", type, type, index);
    FILE* stream = fopen(name, "r");
    fread(out, len, 1, stream);
    fclose(stream);
}

void decode_block(int len, int e, u8** out_a, u8** out_b) {
    int i, j;
    if (e>=K) { 
        u8* stripe_a[K];
        u8* stripe_b[K];
        for (i=0; i<K; i++) {
            stripe_a[i]=malloc(len);
            get_block(len, "a", i, stripe_a[i]);
            stripe_b[i]=malloc(len);
            get_block(len, "b", i, stripe_b[i]);
        }
        hhk_encode_parity(len, e, stripe_a, stripe_b, out_a, out_b);
    } else {
        // 1. decode b
        u8* stripe_b_left[K-1];
        j=0;
        for (i=0; i<K; i++) {
            if (i!=e) {
                stripe_b_left[j]=malloc(len);
                get_block(len, "b", i, stripe_b_left[j++]);
            }
        }
        u8 fb[len];
        get_block(len, "pb", 0, fb);
        hhk_decode_e1(len, e, K, stripe_b_left, fb, out_b);

        // gen full stripe_b
        u8* stripe_b[K];
        for (i=0; i<K; i++) {
            if (i<e) stripe_b[i] = stripe_b_left[i];
            else if (i==e) stripe_b[i] = out_b[0];
            else stripe_b[i] = stripe_b_left[i-1];
        }

        // 2. decode a
        if (e==K-1) {
            u8* parity2=malloc(len), *parity3=malloc(len), *parity4=malloc(len);
            get_block(len, "pb", 2, parity2);
            get_block(len, "pb", 3, parity3);
            get_block(len, "pa", 1, parity4);
            hhk_decode_a_last(len, stripe_b, parity2, parity3, parity4, out_a);
        } else { // decode data 0~K-1
            u8* stripe_a_left[2];
            j=0;
            for (i=e/3*3; i<(e/3+1)*3; i++) { // 0,1,2->(0,3)  3,4,5->(3,6)  6,7,8->(6,9)
                if (i!=e) {
                    stripe_a_left[j]=malloc(len);
                    get_block(len, "a", i, stripe_a_left[j++]);
                }
            }
            u8* parity=malloc(len);
            get_block(len, "pb", e/3+1, parity);
            hhk_decode_a(len, e, stripe_b, stripe_a_left, parity, out_a);
        }
    }
}

int main(int argc, char* argv[])
{
    int e=0;
    if (argc==2) {
        e = atoi(argv[1]);
    }

    encode_matrix = malloc(M * K);
    gf_gen_cauchy1_matrix(encode_matrix, M, K);

    // TODO: 内存泄漏，文件描述符释放

    encode_file("a.txt");
    int len=4;
    printf("----- test -----\n");

    u8* out_a[1];
    u8* out_b[1];
    out_a[0]=malloc(len);
    out_b[0]=malloc(len);

    printf("e: %d\n", e);
    decode_block(len, e, out_a, out_b);
    printf("%d: (%s) (%s)\n", e, out_a[0], out_b[0]);
}
