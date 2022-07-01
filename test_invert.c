#include <stdio.h>
#include <stdlib.h>
#include <isa-l/erasure_code.h>
#include "helper.h"

typedef unsigned char u8;

int main() {
    // 1 0 0 0      0 1 0 0
    // 0 1 0 0      0 0 1 0
    // 0 0 0 1      0 0 0 1
    // 4 6 7 8      4 6 7 8

    /*

        1 0 0 0   b1   b1
        0 1 0 0   b2   b2
     x  0 0 1 0 * b3 = b3 x
        0 0 0 1   b4   b4
        4 6 7 8        f(b)

    out1:
    1, 0, 0, 0,             b1       b1
    0, 1, 0, 0,             b2       b2
    210, 187, 185, 186, *   b4   =   x (应该是 b3) 
    0, 0, 1, 0,             f(b)     b4

    out2:
    143, 200, 2, 71,        b2       x (应该是 b1)
    1, 0, 0, 0,             b3       b2
    0, 1, 0, 0,         *   b4   =   b3
    0, 0, 1, 0,             f(b)     b4
    */

    // 规律：丢失哪一个，生成的系数行就在哪个位置。其他位置都是单位行获取原本的bi
    // 因此只需要这个单独的系数行。

    u8 in1[16] = {1,0,0,0, 0,1,0,0, 0,0,0,1, 4,6,7,8};
    u8 in2[16] = {0,1,0,0, 0,0,1,0, 0,0,0,1, 4,6,7,8};
    u8 out1[16], out2[16];
    show_matrix(in1, 4, 4, "in1");
    gf_invert_matrix(in1, out1, 4);
    show_matrix(out1, 4, 4, "out1");

    show_matrix(in2, 4, 4, "in2");
    gf_invert_matrix(in2, out2, 4);
    show_matrix(out2, 4, 4, "out2");

    printf("----- verify invert matrix recover data of BYTE -----\n");
    // 验证 out1 能否恢复 b3
    u8 b1=33, b2=213, b3=58, b4=9;
    printf("init b3: %d\n", b3);
    u8 fb = gf_mul(4,b1)^gf_mul(6,b2)^gf_mul(7,b3)^gf_mul(8,b4);
    printf("fb: %d\n", fb); // 174
    // 210, 187, 185, 186 * b1 b2 b4 f(b)
    u8 recover_b3 = gf_mul(210,b1)^gf_mul(187,b2)^gf_mul(185,b4)^gf_mul(186,fb);
    printf("recover b3: %d\n", recover_b3); // 58

    printf("----- verify ec_encode_data -----\n");
    u8 encode_matrix[] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1, 4,6,7,8};
    u8 gftbls[1*4*32];
    ec_init_tables(4, 1, &encode_matrix[16], gftbls);
    
    u8 b_data1[] = {33,213,58,9};
    u8 b_data2[] = {73,99,108,232};
    u8* data[4];
    for (int i=0; i<4; i++) { // len=2, 2 bytes / block
        data[i] = malloc(2);
        data[i][0] = b_data1[i];
        data[i][1] = b_data2[i];
    }
    u8* coding[1];
    coding[0] = malloc(2);

    ec_encode_data(2, 4, 1, gftbls, data, coding);
    u8 fb1=0, fb2=0;
    u8 coefficients[] = {4,6,7,8};
    for (int i=0; i<4; i++) {
        fb1 ^= gf_mul(coefficients[i], b_data1[i]);
        fb2 ^= gf_mul(coefficients[i], b_data2[i]);
    }
    printf("fb1: %d\n", fb1); // 174
    printf("fb2: %d\n", fb2); // 100
    printf("coding[0]: (%d,%d)\n", coding[0][0], coding[0][1]); // (174,100) 与普通算法结果相同（向量点积，加是异或，乘是gf_mul）
    // ec_encode_data配合gftbls, 对于块大小len为多个字节，可以并行编码有性能提升。
    // ec_base.c/ec_encode_data_base 其实就是对块里的每个字节，一个个乘以系数生成的gf_tbls, 编码成对应的字节，写到输出块对应位置。
}
