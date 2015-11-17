/*------------------------------------------------------------------*/
/* UTILITY-LIB.C                                                    */
/*                                                                  */
/* Version 0.1, 08.11.95, Juergen Schreuer                          */
/*                                                                  */
/*------------------------------------------------------------------*/

/*-- INCLUDE -------------------------------------------------------*/

#if !defined( _UTILITY_LIB )
#define _UTILITY_LIB
#define _UL_DEBUG_LEVEL       0

#include <stdio.h>    
#include <stdlib.h>   
#include <string.h>   
#include <signal.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <time.h>


/*-- DEFINES -------------------------------------------------------*/

/* Error codes */

#define    FOPEN_ERROR          1001  /* Unable to open file */
#define    FILE_EMPTY           1002  /* File contains no usable data */
#define    FUNC_MODE_UNKNOWN    1003  /* Unknown function mode */

/* Control flags */

#define    M_ON               1  
#define    M_OFF              0

#define    M_TRUE             1
#define    M_FALSE            0

#define    M_SET              1
#define    M_UNSET            0

#define    M_VALID            1
#define    M_INVALID          0

#define    M_SUCCESS          0
#define    M_ERROR            1
#define    M_ABORT            2
#define    M_CONTINUE         3
#define    M_WARNING          4
#define    M_NOTICE           5

/* Constants */

#define    LINE_BUF_SIZE      500 

/* Token types */

#define    TOK_INVALID        0
#define    TOK_USR_COMMENT    1
#define    TOK_C_COMMENT      2
#define    TOK_ALNUM          3
#define    TOK_STRING         4
#define    TOK_CONTROL        5
#define    TOK_EOS            6
  
/*-- STRUCTURE DEFINITIONS -----------------------------------------*/

typedef struct               /* Token structure */
    {
      char *actual_token;    /* Pointer to current token */
      char *next_token;      /* Pointer to next token */
      int   tok_type;        /* Current token type */
      int   tok_length;      /* Length of current token */
    } TOKEN;
    
typedef struct               /* Tokentype structure */
    {
      int usr_comment;       
      int c_comment;
      int alnum;
      int string;
      int control;
      int invalid;
    } TOKTYPES;  

/*-- GLOBAL VARIABLES ----------------------------------------------*/

int   ErrorStatus;                 /* Error status variable */
char  ErrorMessage[LINE_BUF_SIZE]; /* Error message buffer */
int   LineCounter;                 /* Line counter used by fgetline() */
/*FILE *LOGstream = stderr;*/  /* Stream used for writing log information */
/*FILE *INPstream = stdin;*/   /* Stream used for standard read */
/*FILE *OUTstream = stdout;*/  /* Stream used for standard write */

/*-- PROTOTYPES ----------------------------------------------------*/

static char  *char_insert( char *insert_here, char c );
int           chk_cmd( char *cmd, char **cmd_list, int mode );
int           count_tokens( char *s, TOKTYPES *tt ); 
int           error_msg( char *msg, int err_handle );
int           fgetline( FILE *stream, char *str );
char         *get_time( char *format );
unsigned int  getbits( unsigned int x, int p, int n );
unsigned int  maskbits( unsigned int x, unsigned int m, int mode );
unsigned int  setbitOFF( unsigned int x, int n );
unsigned int  setbitON( unsigned int x, int n );
char         *str_token( TOKEN *tok );

/*-- Functions --------------------------------------------------------*/

/*---------------------------------------------------------------------*/
/* unsigned int setbitON( unsigned int x, int n );                     */
/* unsigned int setbitOFF( unsigned int x, int n );                    */
/* usigned int getbits( unsigned int x, int p, int n );                */
/* unsigned int maskbits( unsigned int x, unsigned int m, int mode );  */
/*                                                                     */
/* Bit manipulations:                                                  */
/*   setbitON : Turns on bit n of x.                                   */
/*   setbitOFF: Turns off bit n of x.                                  */
/*   getbits  : Get n bits from position p of x.                       */
/*   maskbits : All bits set to one in mask m are turned on ( M_ON )   */
/*              or off ( M_OFF ) in x.                                 */
/*---------------------------------------------------------------------*/

unsigned int setbitON( unsigned int x, int n )
{
  return( x | ( 1 << n ) );
}

unsigned int setbitOFF( unsigned int x, int n )
{
  return( ( x | ( 1 << n ) ) & ~( 1 << n ) );
}

unsigned int getbits( unsigned int x, int p, int n )
{
  return( ( x >> ( p + 1 - n ) ) & ~( ~0 << n ) );
}

unsigned int maskbits( unsigned int x, unsigned int m, int mode )
{
  if( mode == M_ON ) return( x | m );
  if( mode == M_OFF ) return( ( x | m ) & ~m );
  return( 0 );
}

/*---------------------------------------------------------------------*/
/* int count_tokens( char *s, TOKTYPES *tt );                          */
/*                                                                     */
/* Counts the tokens in string s and returns the total number of tokens*/
/* found. If tt != NULL the structure tt contains the number of found  */
/* tokens of the known token types.                                    */
/*---------------------------------------------------------------------*/

int count_tokens( char *s, TOKTYPES *tt )
{
  int    count = 0;
  char  *scp;
  TOKEN  tok;
  
  scp = (char *)malloc( strlen( s ) * sizeof(char));
  strcpy( scp, s );
  if( NULL != tt )
  {
    tt->usr_comment = 0;
    tt->c_comment = 0;
    tt->alnum = 0;
    tt->string = 0;
    tt-> control = 0;
    tt->invalid = 0;
  }
  tok.actual_token = tok.next_token = scp;
  while( 1 )
  {
    str_token( &tok );
    if( tok.tok_type == TOK_EOS ) break;
    switch( tok.tok_type )
    {
      case TOK_ALNUM:
        if( NULL != tt ) tt->alnum++;
        count++;
        break;
      case TOK_USR_COMMENT:
        if( NULL != tt ) tt->usr_comment++;
        count++;
        break;
      case TOK_C_COMMENT:
        if( NULL != tt ) tt->c_comment++;
        count++;
        break;
      case TOK_STRING:
        if( NULL != tt ) tt->string++;
        count++;
        break;
      case TOK_CONTROL:
        if( NULL != tt ) tt->control++;
        count++;
        break;
      case TOK_INVALID:
        if( NULL != tt ) tt->invalid++;
        count++;
        break;
      default:
        break;
    }
  }
  free( scp );
  return( count );
}

/*----------------------------------------------------------------------------*/
/* char *get_time( char *format );                                            */
/*                                                                            */
/* Gets the current date and time, formats it according to format and returns */
/* a pointer to the resulting string.                                         */
/*----------------------------------------------------------------------------*/


char *get_time( char *format )
{
  static char CurTime[32];
  time_t      t;
  
  if( -1 == time( &t ) || 0 == strftime( CurTime, 31, format, localtime( &t )))
  {
    CurTime[0] = '?';
    CurTime[1] = '\0';
  }
  
  return( CurTime );
}

/*---------------------------------------------------------------------*/
/* int fgetline( FILE *stream, char *str );                            */
/*                                                                     */
/* Reads a maximum of LINE_BUF_SIZE characters from stream until '\n', */
/* '\0' or EOF is reached and returns the number of characters read.   */
/* Until EOF is reached every call to fgetline increments the global   */
/* variable LineCounter                                                */
/*---------------------------------------------------------------------*/

int fgetline( FILE *stream, char *str )
{
  int c, i = 0;
  while( i < LINE_BUF_SIZE && ( c = fgetc( stream ) ) != EOF && c != '\n' && c != '\n' )
  {
    *(str + i++) = (char) c;
  }
  *(str+i) = '\0';
  if( i == 0 && c == EOF )
    return( EOF );
  LineCounter ++;
  return( i );
}

/*----------------------------------------------------------------------*/
/*  char *str_token( TOKEN *tok );                                      */
/*                                                                      */
/*  Searches for tokens in string. Before the first call to str_token   */
/*  tok->actual_token and tok->next_token must be initialized with      */
/*  the pointer to the beginning of string. The entire string can be    */
/*  parsed by subsequent calls to str_token. tok->tok_type contains     */
/*  the token type of the actual token, tok->tok_length its length.     */
/*                                                                      */
/*  Blanks ' ' and tabs '\t' are interpreted as token delimiters.       */
/*  If a token is found, '\0' is inserted in the string at end of token,*/
/*  the internal pointer next_token is set to the next character        */
/*  following the actual token and a pointer to the actual token is     */
/*  returned.                                                           */
/*                                                                      */
/*  Token types:                                                        */
/*    TOK_INITIALIZED  : Internal pointers are successfully initialized,*/
/*                       no token search is performed.                  */
/*                       Return value: pointer to string                */
/*    TOK_USR_COMMENT  : '#' is first symbol of token, rest of string   */
/*                       is interpreted as one token.                   */
/*                       Return value: pointer to token.                */
/*    TOK_C_COMMENT    : '/','*' are the first two characters of the    */
/*                       token, which ends at the following '*','/'     */
/*                       sequence or the end of string.                 */
/*                       Return value: pointer to token.                */
/*    TOK_ALNUM        : Alpha-numerical token, consisting of printable */
/*                       characters ('\x21' - '\x7E').                  */
/*                       Return value: pointer to token.                */
/*    TOK_STRING       : '"' is first character of token, which ends at */
/*                       the next '"' or the end of string. The second  */
/*                       '"' is replaced by '\0'.                       */
/*                       Return value: pointer to second character of   */
/*                       token, i.e. character following the first '"'. */ 
/*    TOK_INVALID      : actual character < '\x20' or > '\x7E', except  */
/*                       '\0' and '\t', no token is defined!            */
/*                       Return value: pointer to character, no '\0' is */
/*                       inserted!                                      */
/*    TOK_EOS          : End Of String is reached, no token is defined. */
/*                       Return value: NULL.                            */
/*----------------------------------------------------------------------*/

char *str_token( TOKEN *tok )
{  
  /* skip leading blanks and tabs */

  while( *(tok->next_token) == ' ' || *(tok->next_token) == '\t' ) (tok->next_token)++;  
  tok->actual_token = tok->next_token;

  /* search for open c-comment */

  while( *++(tok->next_token) )
    if( *(tok->next_token) == '/' )
    {
      if( *(tok->next_token-1) == '*' )
      {
        tok->tok_type = TOK_C_COMMENT;
        tok->next_token = char_insert( tok->next_token + 1, '\0' );
        tok->tok_length = (int)(tok->next_token - tok->actual_token - 1);
        return( tok->actual_token );
      }
      if( *(tok->next_token+1) == '*' )
        break;
    }    

  tok->next_token = tok->actual_token;    
  
  /* *next_token is printable character, except space */

  if( '\x20' < *(tok->next_token) && '\x7F' > *(tok->next_token) )
    switch( *(tok->next_token) )
    {
      case '#' :                        /* token is user comment */
        tok->tok_type = TOK_USR_COMMENT;
        while( *++(tok->next_token) );
        tok->tok_length = (int)(tok->next_token - tok->actual_token);
        return( tok->actual_token );
      case '"' :                        /* token is string, delimiter '"' */
        tok->tok_type = TOK_STRING;
        while( *++(tok->next_token) )
          if( *(tok->next_token) == '"' )
          {
            *(tok->next_token)++ = '\0';
            break;
          }
        (tok->actual_token)++;
        tok->tok_length = (int)(tok->next_token - tok->actual_token);
        if( *(tok->next_token) ) (tok->tok_length)--;
        return( tok->actual_token );       
      case '\'' :                       /* token is string, delimiter '\'' */
        tok->tok_type = TOK_STRING;
        while( *++(tok->next_token) )
          if( *(tok->next_token) == '\'' )
          {
            *(tok->next_token)++ = '\0';
            break;
          }
        (tok->actual_token)++;
        tok->tok_length = (int)(tok->next_token - tok->actual_token);
        if( *(tok->next_token) ) (tok->tok_length)--;
        return( tok->actual_token );       
      case '/' :
        if( *++(tok->next_token) == '*' )      /* token is c comment */
        {
          tok->tok_type = TOK_C_COMMENT;
          while( *++(tok->next_token) )
            if( *(tok->next_token) == '*' && *(tok->next_token+1) == '/' )
            {
              tok->next_token = char_insert( tok->next_token + 2, '\0' );
              break;
            }  
          tok->tok_length = (int)(tok->next_token - tok->actual_token);
          if( *(tok->next_token) ) (tok->tok_length)--;
          return( tok->actual_token );
        }
        (tok->next_token)--;
      default :                         /* token is alphanumerical */    
        tok->tok_type = TOK_ALNUM;
        while( isgraph( *++(tok->next_token) ) );
        tok->next_token = char_insert( tok->next_token, '\0' );
        tok->tok_length = (int)(tok->next_token - tok->actual_token - 1);
        return( tok->actual_token );
      }
  else                                    /* token is not printable */
    switch( *(tok->next_token) )
    {
      case '\0' :
        tok->tok_type = TOK_EOS;          /* End Of String reached */
        return( NULL );
      default   :
        if( *(tok->next_token) <= '\x20' )
          tok->tok_type = TOK_CONTROL;   /* Control character */
        else
          tok->tok_type = TOK_INVALID;   /* Invalid token */
        (tok->next_token)++;
        tok->tok_length = 1;
        return( NULL );   
     }
}

/*----------------------------------------------------------------------*/
/* static char *char_insert( char *insert_here, char c );               */
/*                                                                      */
/* Inserts char c at pointer insert_here in a string terminated         */
/* with '\0' and returns pointer to next character in string.           */
/*----------------------------------------------------------------------*/

static char *char_insert( char *insert_here, char c )
{
  char *ptr = insert_here;

  while( *ptr++ );
  while( insert_here < ptr-- ) *(ptr+1) = *ptr;
  *insert_here = c;
  return( insert_here + 1 );
}  

/*-----------------------------------------------------------------------*/
/* int chk_cmd( char *cmd, char **cmd_list, int mode );                  */
/*                                                                       */
/* Checks, if cmd is a command included in the list cmd_list. Returns    */
/* the number (>0) of the command or 0 if it is unknown.                 */
/*                                                                       */
/* mode == 'u' : convert all characters of cmd to uppercase.             */
/*      == 'l' : convert all characters of cmd to lowercase.             */
/*-----------------------------------------------------------------------*/

int chk_cmd( char *cmd, char **cmd_list, int mode )
{
  int   ncmd;
  char *ptr;

  ptr  = cmd - 1;
  switch( mode )
  {
    case 'u' :
      while( *++ptr )
        *ptr = toupper( *ptr );
      break;
    case 'l' :  
      while( *++ptr )
        *ptr = tolower( *ptr );
      break;
    default:
      break;
  }

  ncmd = atoi( *cmd_list );
  while( ncmd > 0 && strcmp( cmd, cmd_list[ncmd] ) ) ncmd--;

  return( ncmd );
}   

/*----------------------------------------------------------------------*/
/* int error_msg( char *msg, int err_handle );                          */
/*                                                                      */
/* Simple error handler, which writes a message msg to LOGstream.       */
/* err_handle == M_ABORT     : Write error message and force normal     */
/*                             program termination with return value    */
/*                             EXIT_FAILURE.                            */
/*               M_CONTINUE  : Write error message and continue program */
/*                             execution.                               */
/*               M_WARNING   : Write warning message and continue.      */
/*               M_NOTICE    : Write notice and continue.               */
/*----------------------------------------------------------------------*/
/*
int error_msg( char *msg, int err_handle )
{
  switch( err_handle )
  {
    case M_ABORT :
      fprintf( LOGstream, "FATAL ERROR - %s\nprogram termination\n", msg );
      exit( EXIT_FAILURE );
      break;
    case M_CONTINUE :
      fprintf( LOGstream, "ERROR - %s\n", msg );
      return( M_ERROR );
    case M_WARNING :
      fprintf( LOGstream, "WARNING - %s\n", msg );
      return( M_WARNING );
    case M_NOTICE :
      fprintf( LOGstream, "NOTICE - %s\n", msg );
      return( M_NOTICE );
    default :
      break;
  }    
  return( M_ERROR );
}     
*/
/*----- End of library code --------------------------------------------*/

/*======================================================================*/

/*-- Debug main --------------------------------------------------------*/

#if _UL_DEBUG_LEVEL > 0 

int main( int argc, char **argv )
{
  int   tt, n = 0;
  char *ptr;
  TOKEN testtok;
  
  printf( "\n\nUTILITY_LIB - DEBUG MODE\n\n" );

  if( argc != 2 )
  {
    printf( "usage: utility_lib \"test string\"\n" );
    return( -1 );
  }
    
  printf( "Test string: <%s>\n\n", *(argv+1) );
  testtok.actual_token = testtok.next_token = *(argv+1);
  while( NULL != ( ptr = str_token( &testtok )))
    printf( "#%d : length=%4d token_type=%2d token=<%s>\n", ++n, testtok.tok_length, testtok.tok_type, ptr );

  return( 0 );
}

#endif
#endif

