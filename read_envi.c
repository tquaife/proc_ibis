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
    int d_type ;
    char interleave[5];
    float wavels[5000];
} envi_meta ;

void read_envi_header(FILE * ifp, envi_meta * hdr)
{

    char line[MAX_LINE_LEN] ;
    char tmpc1[MAX_LINE_LEN] ;
    char tmpc2[MAX_LINE_LEN] ;
    char *tmpsp;
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
        else if( !strncmp( tmpc1, "data type", 9) ) hdr->d_type = atoi( tmpc2 ); 
        
        /*wavelength is a special case*/
        else if( !strcmp( tmpc1, "Wavelength ")){
            while( fgets( line, MAX_LINE_LEN, ifp ) ){
                tmpsp = strtok(line, ",");
                while( tmpsp != NULL ) {
                     /*check for the end of the wavelength record*/
                     if(strchr(tmpsp, '}')!=NULL){break_out=TRUE;break;}
                     if(!is_string_blank(tmpsp))
                         /*if we get here we should have a valid wavelength*/
                         *(hdr->wavels+nwvls++) = atof(tmpsp);
                     /*get next token*/
                     tmpsp = strtok(NULL, ",");                         
                 }
                /*if break_out==TRUE we are at the end of the wavelength section*/     
                if( break_out ) break ;
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
    for(int i=0; i<hdr.n_bands; i++) printf("%0.2f ", hdr.wavels[i]);
    printf("\n");
    printf("samples: %d\n",hdr.n_samples);
    printf("lines: %d\n",hdr.n_lines);
    printf("bands: %d\n",hdr.n_bands);
    printf("interleave: %s\n",hdr.interleave);
    printf("data type: %d\n",hdr.d_type);
 
    return(EXIT_SUCCESS);
}
