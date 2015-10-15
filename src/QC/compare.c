#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

void write_index(FILE* file, unsigned long start, unsigned long end, char* chromo_name);

int main(int argc, char **argv) {

	if (argc != 3) {
		printf("Usage: <file1> <file2>\n");
		exit(EXIT_FAILURE);	
	}

	FILE *file1, 
		 *file2,
		 *only_f1,
		 *only_f2,
		 *overlaps;

	if ((file1 = fopen(argv[1], "r")) == NULL) {
		printf("ERROR: Could not open first file\n");
		exit(EXIT_FAILURE);
	}

	if ((file2 = fopen(argv[2], "r")) == NULL) {
		printf("ERROR: Could not open second file\n");
		exit(EXIT_FAILURE);
	}

	size_t size_f1_name = strlen(argv[1]);
	size_t size_f2_name = strlen(argv[2]);
	char name_cat[(size_f1_name + size_f2_name + 6)];
	strcpy(name_cat, argv[1]);
	strcat(name_cat, "_x_");
	strcat(name_cat, argv[2]);
		
	char overlaps_name[(size_f1_name + size_f2_name + 15)];
	strcpy(overlaps_name, name_cat);
	strcat(overlaps_name, "_overlaps");
	if ((overlaps = fopen(overlaps_name, "w")) == NULL) {
		printf("Could not create overlaps file\n");
		exit(EXIT_FAILURE);
	}
	
	char only_f1_name[(size_f1_name + size_f2_name + 14)];
	strcpy(only_f1_name, name_cat);
	strcat(only_f1_name, "_only_f1");
	if ((only_f1 = fopen(only_f1_name, "w")) == NULL) {
		printf("Could not create only_f1 file\n");
		exit(EXIT_FAILURE);
	}

	char only_f2_name[(size_f1_name + size_f2_name + 14)];
	strcpy(only_f2_name, name_cat);
	strcat(only_f2_name, "_only_f2");
	if ((only_f2 = fopen(only_f2_name, "w")) == NULL) {
		printf("Could not create only_f2 file\n");
		exit(EXIT_FAILURE);
	}

	unsigned long other = 0, 
		x_in_file1 = 0, 
		x_in_file2 = 0, 
		x_both = 0;

	char c, 
		 d;

	/* Comparing base by base */

	while ((c = getc(file1)) != EOF) {
		if ((d = getc(file2)) == EOF) {
			printf("ERROR: Input files different lengths\n");
			exit(EXIT_FAILURE);
		}

		if (c == '>') {
			if (d != '>') {
				printf("ERROR: Discrepancy in chromosome lengths\n");
				exit(EXIT_FAILURE);
			}
			else {
				while ((c = getc(file1)) != '\n') {
					d = getc(file2);
					if (c == EOF || d == EOF) {
						printf("ERROR: Unexpected end of file\n");
						exit(EXIT_FAILURE);
					}
					if (d == '\n') {
						printf("ERROR: Discrepancy in chromosome names\n");
						exit(EXIT_FAILURE);
					}
				}

				if ((d = getc(file2)) == EOF) {
					printf("ERROR: Unexpected end of file in file2\n");
					exit(EXIT_FAILURE);
				}

				if (d != '\n') {
					printf("ERROR: Discrepancy in chromosome names\n");
					exit(EXIT_FAILURE);
				}

				if (((c = getc(file1)) == EOF) || ((d = getc(file2)) == EOF)) {
					printf("ERROR: Unexpected end of file\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		if (c != '\n') {
			if (d == '\n') {
				printf("c = %c\n", c);
				printf("ERROR: Unexpected newline in file2\n");
				exit(EXIT_FAILURE);
			}
			else if (c == 'X' && d == 'X') 
				x_both++;
			else if (c == 'X' && d != 'X') 
				x_in_file1++;
			else if (c != 'X' && d == 'X') 
				x_in_file2++;
			else if (c != 'X' && d != 'X') {
				other++;
			}
		}

		else if (c == '\n' && d != '\n') {
			printf("d = %c\n", d);
			printf("ERROR: Unexpected newline in file1\n");
			exit(EXIT_FAILURE);
		}
	}
	
	unsigned long total = x_both + x_in_file1 + x_in_file2 + other;
	printf("Agree: %lu\nIn f1 but not f2: %lu\nIn f2 but not f1: %lu\n\nTotal bases: %lu\n\n", x_both, x_in_file1, x_in_file2, total);



	/* Check to see where there is overlap between shredded reads and peak k-mers */

	rewind(file1);
	rewind(file2);

	int state; /* 0: neither Xs; 1: f1 X but not f2; 2: f2 X but not f1; 3: both Xs */
	unsigned long count_overlap = 0,
		count_f1 = 0,
		count_f2 = 0,
		base_index = 1,
		start;

	char* chromo_name;
	if ((chromo_name = malloc(60 * sizeof(char))) == NULL) {
		printf("Out of memory\n");
		exit(EXIT_FAILURE);
	}

	while ((c = getc(file1)) != EOF) {
		if ((d = getc(file2)) == EOF) {
			printf("ERROR: Input files different lengths\n");
			exit(EXIT_FAILURE);
		}

		if (c != '\n' && c != '>') {
			if (state == 0) {
				if (c == 'X' && d != 'X') {
					state = 1;
					count_f1++;
					start = base_index;
				}
				else if (c != 'X' && d == 'X') {
					state = 2;
					count_f2++;
					start = base_index;
				}
				else if (c == 'X' && d == 'X') {
					state = 3;
					count_overlap++;
					start = base_index;
				}
			} 
			
			else if (state == 1) {
				if (c == 'X' && d == 'X') {
					state = 3;
					count_overlap++;
					write_index(only_f1, start, base_index, chromo_name);
					start = base_index;
				}
				else if (c != 'X' && d == 'X') {
					state = 2;
					count_f2++;
					write_index(only_f1, start, base_index, chromo_name);
					start = base_index;
				}
				else if (c != 'X' && d != 'X') {
					state = 0;
					write_index(only_f1, start, base_index, chromo_name);
				}
			}

			else if (state == 2) {
				if (c == 'X' && d == 'X') {
					state = 3;
					count_overlap++;
					write_index(only_f2, start, base_index, chromo_name);
					start = base_index;
				}
				else if (c == 'X' && d != 'X') {
					state = 1;
					count_f1++;
					write_index(only_f2, start, base_index, chromo_name);
					start = base_index;
				}
				else if (c != 'X' && d != 'X') {
					state = 0;
					write_index(only_f2, start, base_index, chromo_name);
				}
			}

			else if (state == 3) {
				if (c != 'X' && d != 'X') {
					state = 0;
					write_index(overlaps, start, base_index, chromo_name);
				}
				else if (c == 'X' && d != 'X') {
					state = 1;
					count_f1++;
					write_index(overlaps, start, base_index, chromo_name);
					start = base_index;
				}
				else if (c != 'X' && d == 'X') {
					state = 2;
					count_f2++;
					write_index(overlaps, start, base_index, chromo_name);
					start = base_index;
				}
			}
			base_index++;
		}

		/* Start of new chromosome => Get name */
		else if (c == '>') {
			free(chromo_name);
			int chromo_name_buffsize = 60;
			if ((chromo_name = malloc(chromo_name_buffsize * sizeof(char))) == NULL) {
				printf("Out of memory\n");
				exit(EXIT_FAILURE);
			}
			int chromo_name_index = 0;
			while ((c = getc(file1)) != '\n') {
				if ((d = getc(file2)) != c || c == EOF) {
					printf("ERROR: Discrepancy in chromosome names between files (files should be identical)\n");
					exit(EXIT_FAILURE);
				}
				if (chromo_name_index == (chromo_name_buffsize - 1)) {
					chromo_name_buffsize *= 2;
					char* tmp = realloc(chromo_name, chromo_name_buffsize * sizeof(char));
					if (tmp != NULL) {
						chromo_name = tmp;
					}
					else {
						printf("ERROR: Out of memory\n");
						free(chromo_name);
						exit(EXIT_FAILURE);
					}
				}
				chromo_name[chromo_name_index++] = c;
			}

			if ((d = getc(file2)) == EOF) {
				printf("ERROR: File2 too short\n");
			}
			chromo_name[chromo_name_index] = '\0';
			base_index = 1;

			/* Find correct state for start of new chromosome */

			if (c != 'X' && d != 'X') {
				state = 0;
			}
			else if (c == 'X' && d != 'X') {
				state = 1;
			}
			else if (c != 'X' && d == 'X') {
				state = 2;
			}
			else if (c == 'X' && d == 'X') {
				state = 3;
			}
			else {
				printf("File incorrectly formatted\n");
				exit(EXIT_FAILURE);
			}

		}
	}
	if (state == 1) {
		write_index(only_f1, start, base_index, chromo_name);
	}
	else if (state == 2) {
		write_index(only_f2, start, base_index, chromo_name);
	}
	else if (state == 3) {
		write_index(overlaps, start, base_index, chromo_name);
	}

	printf("overlap: %lu\nf1: %lu\nf2: %lu\n\n", count_overlap, count_f1, count_f2);	

	free(chromo_name);
	fclose(file1);
	fclose(file2);
	fclose(only_f1);
	fclose(only_f2);
	fclose(overlaps);

	return 0;
}


void write_index(FILE* file, unsigned long start, unsigned long end, char* chromo_name) {
	size_t total_size = strlen(chromo_name) + 60;
	char out_str[total_size];
	snprintf(out_str, total_size, "%-20s %12lu %12lu %12lu\n", chromo_name, start, (end - 1), (end - start));
	fputs(out_str, file);

	return;
}	

