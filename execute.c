#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "field.h"
#include "functions.h"
  
void main(void)
{
	int option=4;
	char input[100], encoded[100], decoded[100];

	printf("\n");
	printf("\tReed-Solomon Codec - GF(256)\n");
	printf("\t----------------------------\n");

	//printf("\n>> Enter the message length: ");
	//scanf("%d",&MSG);
	MSG = 9;
	POWER = 4;						// Galois Field Power
	FIELD = (unsigned int)pow(2,POWER);
	PAR = FIELD-MSG-1;				// 2t: t-error correcting

	// Generate the Zechs Logarithm Tables
	ZechsLogarithms();
	RSDecoder("input.txt","try.txt");
	
/*	printf("\n\t1. Encode\n\t2. Decode\n\t3. Exit\n\n");
	printf(">> Enter your choice: ");
	scanf(" %d",&option);

	while(option!=3)
	{
		switch(option)
		{
		case 1:
			printf("\nEnter the input file: ");
			scanf("%s",input);
			sprintf(encoded,"enc-%s",input);
			RSEncoder(input,encoded);
			option = 4;
			break;
		case 2:
			printf("\nEnter the encoded file: ");
			scanf("%s",encoded);
			printf("\nEnter the decoded filename: ");
			scanf("%s",decoded);
			RSDecoder(encoded,decoded);
			option = 4;
			break;
		case 3:
			printf("\n");
			exit(0);
		default:
			printf("\n\t1. Encode\n\t2. Decode\n\t3. Exit\n\n");
			printf(">> Enter your choice: ");
			fflush(stdin);
			scanf(" %d",&option);
			break;
		};
		
	}*/
	printf("\n");
}

