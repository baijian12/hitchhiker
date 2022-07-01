#ifndef _HELPER_H_
#define _HELPER_H_

/**
 * @brief Display the matrix of size rows*cols, start with "message"
 */
void show_matrix(unsigned char* matrix, int rows, int cols, char* message);

/**
 * @brief Display fragments of size rows*cols, start with "message"
 */
void show_frags(unsigned char** frags, int rows, char* message);

#endif //_HELPER_H_