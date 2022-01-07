#include <stdio.h> /* for printf and scanf */

int main(void)     
{
    char *s, string[101];

    printf("Input a string: ");
    scanf("%100s",string);

    /* loop over characters */
    for (s=string; *s != '\0'; s++) 
	printf("%c\n",*s);
}
