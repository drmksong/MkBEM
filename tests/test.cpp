#include <stdio.h>
#include <string.h>

char *trim_left(char *str);

int main()
{
  char str[256];
  memset(str, '\0', 255);
  strcpy(str, "1  	2	    3	      	4");
  char *str2;
  str2 = str;
  int i = 0, j;
  while (str2)
  {
    sscanf(str2, "%d", &j);
    printf("%d: single=%d \n", i, j);
    printf("all=%s \n", str2);
    str2 = strchr(str2, ' ');
    if (!str2)
      break;
    str2 = trim_left(str2);
    printf("after copy=%s \n", str2);
    i++;
  }
  return 0;
}

char *trim_left(char *str)
{
  while (*str == ' ' || *str == '\t')
    str++;
  return str;
}