#ifndef _HITCHHIKER_H_
#define _HITCHHIKER_H_

#define K 10    // src block numbers
#define P 4     // parity block numbers
#define M 14    // src+parity

typedef unsigned char u8;

/**
 * @brief Fill 0 to a proper size, split to K strpie_a and K stripe_b
 * 
 * @param stream    stream of file
 * @param len       number of bytes in a block. Must be power of 2
 * @param stripe_a  array of pointers of stripe a, size K
 * @param stripe_b  array of pointers of stripe b, size K
 * @returns 0 success, else fail
 */
int split_file(FILE* stream, int len, u8** stripe_a, u8** stripe_b);

/**
 * @brief Encode data using hitchhiker algo
 * 
 * @param len       block size of SUB-STRIPE (half of normal block size)
 * @param stripe_a  input a, size K*len
 * @param stripe_b  input b, size K*len
 * @param parity_a  output parity a, size P*len
 * @param parity_b  output parity b, size P*len
 */
void hhk_encode_data(int len, u8** stripe_a, u8** stripe_b, u8** parity_a, u8** parity_b);

/**
 * @brief Decode substripe 1 error. K-1 data + 1 parity -> data_i
 * 
 * @param len               block size of SUB-STRIPE 
 * @param e                 error data block index (start from 0)
 * @param p                 parity coefficients index in encode_matrix
 * @param substripe_left    K-1 left data, MUST BE ORDERED
 * @param parity            parity block
 * @param out               output buffer
 * @returns 0 successful
 */
int hhk_decode_e1(int len, int e, int p, u8** substripe_left, u8* parity, u8** out);

/**
 * @brief Encode one parity block from substripe
 * 
 * @param len           block size of SUB_STRIPE
 * @param p             parity coefficients index in encode_matrix
 * @param substripe     input data
 * @param out           output buffer
 */
void hhk_encode_1(int len, int p, u8** substripe, u8** out);

/**
 * @brief Decode stripe_a of 1 error
 * 
 * @param len               block size of SUB_STRIPE
 * @param e                 error data block index (start from 0)
 * @param stripe_b          K stripe b
 * @param stripe_a_left     stripe a left in same group 
 * @param parity            corresponding parity
 * @param out               output buffer
 */
void hhk_decode_a(int len, int e, u8** stripe_b, u8** stripe_a_left, u8* parity, u8** out);

/**
 * @brief Decode last data block
 * 
 * @param len       block size of SUB_STRIPE
 * @param stripe_b  K stripe b
 * @param parity2   4,5,6
 * @param parity3   7,8,9
 * @param parity4   4,...,10
 * @param out       output buffer
 */
void hhk_decode_a_last(int len, u8** stripe_b, u8* parity2, u8* parity3, u8* parity4, u8** out);

/**
 * @brief Encode 1 parity of hhk algo
 * 
 * @param len       block size of SUB_STRIPE
 * @param p         parity coefficients index in encode_matrix
 * @param stripe_a  K a input
 * @param stripe_b  K b input
 * @param out_a     out parity a buffer 
 * @param out_b     out parity b buffer
 */
void hhk_encode_parity(int len, int p, u8** stripe_a, u8** stripe_b, u8** out_a, u8** out_b);


#endif //_HITCHHIKER_H_