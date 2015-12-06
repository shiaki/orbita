
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


/* printing functions */

int
orbita_err(char * Msg, ...)
{
  va_list ap;
  va_start(ap, Msg);

  char * line_header = "\033[33m[orbita]\033[31m Error: \033[0m";

  char msg_buffer[256];
  vsprintf(msg_buffer, Msg, ap);

  fprintf(stderr, "%s%s\n", line_header, msg_buffer);

  va_end(ap);

  return 0;
}

int
orbita_warn(char * Msg, ...)
{
  va_list ap;
  va_start(ap, Msg);

  char * line_header = "\033[33m[orbita]\033[32m Warning: \033[0m";

  char msg_buffer[256];
  vsprintf(msg_buffer, Msg, ap);

  fprintf(stderr, "%s%s\n", line_header, msg_buffer);

  va_end(ap);

  return 0;
}

int
orbita_msg(char * Msg, ...)
{
  va_list ap;
  va_start(ap, Msg);

  char * line_header = "\033[33m[orbita]\033[0m ";

  char msg_buffer[256];
  vsprintf(msg_buffer, Msg, ap);

  fprintf(stdout, "%s%s\n", line_header, msg_buffer);

  va_end(ap);

  return 0;
}

int
orbita_dbg(char * Msg, ...)
{
  va_list ap;
  va_start(ap, Msg);

  char * line_header = "\033[33m[orbita]\033[36m Debug: \033[0m";

  char msg_buffer[256];
  vsprintf(msg_buffer, Msg, ap);

  fprintf(stdout, "%s%s\n", line_header, msg_buffer);

  va_end(ap);

  return 0;
}

