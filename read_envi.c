#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "strparse.h"


#define MAX_LINE_LEN 6000

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif


typedef struct  {
    int n_samples ;
    int n_lines ;
    int n_bands ;
    char interleave[5];
    int wavels[5000];
} envi_meta ;


void read_envi_header(FILE * ifp, envi_meta * hdr)
{

    char line[MAX_LINE_LEN] ;
    char tmpc1[MAX_LINE_LEN] ;
    char tmpc2[MAX_LINE_LEN] ;
    int  nlines = 0;
    int  nwvls = 0;
    char break_out=FALSE ;
    size_t pos;
    
    while( fgets( line, MAX_LINE_LEN, ifp ) != NULL ){
        if( is_string_blank( line ) ) continue ;
        nlines++;
        
        /*check it's an ENVI header file*/
        if (nlines==1){
            get_first_string_element( line, tmpc1 );
            if (strcmp(tmpc1,"ENVI")){
                fprintf(stderr,"not an ENVI header file\n");
                exit(EXIT_FAILURE);
            }           
            continue;
        }
        
        /*can only get here if nlines>1 and we're reading
        a valid ENVI header file*/
        
        /*split the line at "=" and tidy up a bit*/
        pos=strcspn(line, "=");
        strncpy(tmpc1, line, pos);
        terminate_string_at_pos(tmpc1,pos);
        copy_end_of_string(tmpc2, line, pos+1 );
        terminate_string_at_newline(tmpc2);
        
        /*get basic ENVI header properties*/
        if( !strncmp( tmpc1, "samples", 7 ) ) hdr->n_samples = atoi( tmpc2 ); 
        else if( !strncmp( tmpc1, "lines", 5 ) ) hdr->n_lines = atoi( tmpc2 ); 
        else if( !strncmp( tmpc1, "bands", 5 ) ) hdr->n_bands = atoi( tmpc2 ); 
        else if( !strncmp( tmpc1, "interleave", 10 ) ) strcpy( hdr->interleave, tmpc2 ); 
        
        /*wavelength is a special case*/
        else if( !strncmp( tmpc1, "Wavelength", 10 ) ){
            pos=strcspn(tmpc2, "{");
            strcpy(line, tmpc2+pos+1);
            //printf("%s",line);
            //printf("%s",tmpc2);
            //exit(-1);
            while( TRUE ){
                while( get_first_string_element( line, tmpc1 ) ){
                    
                    pos=strspn(line, "0123456789.");
                    printf("%d ",pos);
                    /*check for the end of the wavelength record*/
                    pos=strcspn(tmpc1, "}");
                    if(*(tmpc1+pos) == '}'){ break_out=TRUE;break;}
                }
                if( break_out ) break ;
                if (fgets( line, MAX_LINE_LEN, ifp ) == NULL ) break ;
            }

        }        
    }            
    return;
}

int main( int argc, char **argv )
{
 
    FILE *ifp;
    envi_meta hdr;

    if((ifp = fopen( argv[argc-1], "r")) == NULL){
        fprintf( stderr, "%s: unable to open file %s\n", *argv, argv[argc-1]);
        exit(EXIT_FAILURE);
    }

   
    read_envi_header(ifp, &hdr);
    printf("samples: %d\n",hdr.n_samples);
    printf("lines: %d\n",hdr.n_lines);
    printf("bands: %d\n",hdr.n_bands);
    printf("interleave: %s\n",hdr.interleave);
    return(EXIT_SUCCESS);
}
