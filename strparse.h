#ifndef STRPARSE_H
#define STRPARSE_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>

/*
not required:
*//*
#define STRPARSE_MAX_LINE_LEN 50000
*/

char is_string_blank( char * );
char strip_comments( char *, char );
char get_first_string_element( char *, char * );
void terminate_string_at_pos( char *, int );
void copy_end_of_string( char *, char *, int );
void terminate_string_at_newline( char * );
#endif


