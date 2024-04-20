#include<strparse.h>

/*Strip comments from string*/

char	strip_comments( char *l, char c ){

	while( *l != '\0' ){
		if( *l == c ){
		 	*l = '\0' ;
			return 0 ;
		}
		l++ ;	
	}	
	return 1 ;


}


/*add null character at specified point*/
void	terminate_string_at_pos( char *l, int pos ){
        *(l+pos)='\0';
	return;
}

/*add null character at first occurence of c*/
void	terminate_string_at_newline( char *l ){
	while( *l != '\0' ){
		if( *l == '\n' ){
		 	*l = '\0' ;
			return;
		}
		l++ ;	
	}	
	return;
}


/*copy from pos to NULL from src into dest*/
void	copy_end_of_string( char *dest, char *src, int pos ){
        src+=(char) pos;
        while(*src!='\0' ){
           *dest=*src;
           dest++;
           src++;
        }
        *dest='\0';
	return;
}


/*Find if string is empty*/
/*Read K&R section 5.5*/

char is_string_blank( char *s )
{
	
	while( *s != '\0' )
		if( ! isspace( *s++ ) ) 
			return 0 ;
				
	return 1 ;

}


/*
** ------------------------
** get_first_string_element
** ------------------------
*//*
Gets the first group of non whitespace characters,
puts them into element and removes them from the
original string.

version 1.1 bug fixed by adding null termination 
see (*) below.
*/

char	get_first_string_element( line, element )
char	*line;
char	*element;
{

	int	i=0, j=0;

	/* Check that we are not already */
	/* at the end of the line        */

	if( *line == '\0' ) return( 0 );
	
	/* Skip leading white space and */
	/* return zero if end of string */
	/* is encountered               */
	
	while( isspace( *(line+i) ) ) 
		if( *(line+( ++i )) == '\0' ) return( 0 );
	
	
	/* Put the first block of non white 
	   space characters into element */
	
	while( ! isspace( *(line+i) ) ){
		if( ( element[j++] = *(line+(i++) ) ) == '\0' ){
		
			*line = '\0';
			return( 1 );
		
		}
	}
	
	
	/*put a NULL terminator in at j - (*) */
	
	element[ j ] = '\0' ;
	
	/* Copy remainder of line into the beggining */
	/* of itself */
	
	for( j=0; ( *(line+j)=*(line+i) ) != '\0' ; ++j, ++i );
	
	
	
	return( 1 );

}





