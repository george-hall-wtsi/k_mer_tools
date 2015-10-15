#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef struct {
	char* name;
	unsigned long start,
				  end;
} location;

location get_next_loc(FILE* mask_locs);
char* get_chromo_name(FILE* ref);
	
int main(int argc, char** argv) {
	FILE *ref,
		 *mask_locs;

	if (argc != 3) {
		printf("Usage: <reference> <location file>\n");
		exit(EXIT_FAILURE);
	}

	if ((ref = fopen(argv[1], "r")) == NULL) {
		printf("ERROR: Could not open reference file\n");
		exit(EXIT_FAILURE);
	}

	if ((mask_locs = fopen(argv[2], "r")) == NULL) {
		printf("ERROR: Could not open location file\n");
		exit(EXIT_FAILURE);
	}


	int num_locations = 0;	
	char c, /* Counts along masked reference */
		 d; /* Counts along unmasked reference */
	while ((c = getc(mask_locs)) != EOF) {
		if (c == '\n') {
			num_locations++;
		}
	}

	rewind(mask_locs);
	char* chromo_name;
	getc(ref); /* Move pointer onto start of first chromo name, instead of '>' character */
	chromo_name = get_chromo_name(ref);
	if ((d = getc(ref)) == EOF) {
		printf("ERROR: Reference file too short\n");
		exit(EXIT_FAILURE);
	} 
	int iCount = 0;
	unsigned long base_index = 1;

	/* For each location caught in mask */
	while (iCount < num_locations) {
		location loc;
		loc = get_next_loc(mask_locs);
		if (loc.end < loc.start) {
			printf("ERROR: End point of location cannot be less than start point\n");
			exit(EXIT_FAILURE);
		}

		/* Make sure that we are in the correct chromosome */
		if (strcmp(chromo_name, loc.name) != 0)	{
			while (strcmp(chromo_name, loc.name) != 0) {
				while ((d = getc(ref)) != '>');
				free(chromo_name);	
				chromo_name = get_chromo_name(ref);
			}
			base_index = 1;
		}

		/* Navigate to correct location on chromosome */
		while (base_index < loc.start) {
			if ((d = getc(ref)) == EOF) {
				printf("ERROR: Reference file too short [2]\n");
				exit(EXIT_FAILURE);
			}
			else if (d == 'A' || d == 'C' || d == 'G' || d == 'T' || d == 'N') {
				base_index++;
			}
			else if (d == '>') {
				printf("ERROR: Reached end of chromosome without finding location\n");
				exit(EXIT_FAILURE);
			}
		}
			
		size_t str_len_required = ((loc.end - loc.start) + 2);
		char sequence[str_len_required];
		int i = 0;
		while (base_index <= loc.end) {
			if (d == 'A' || d == 'C' || d == 'G' || d == 'T' || d == 'N') {
				sequence[i++] = d;
				base_index++;
			}
			else if (d == '>') {
				printf("ERROR: Reached end of chromosome without finding location\n");
				exit(EXIT_FAILURE);
			}
			if ((d = getc(ref)) == EOF) {
				printf("ERROR: Reference file too short [3]\n");
				exit(EXIT_FAILURE);
			}
		}
		sequence[i] = '\0';
		if ((str_len_required - 1) >= 1) {
			printf(">%s_%lu_%lu\n%s\n", loc.name, loc.start, loc.end, sequence);
		}
		iCount++;
	}
	
	free(chromo_name);
	fclose(ref);
	fclose(mask_locs);
	return 0;
}

char* get_chromo_name(FILE* ref) {
	int chromo_name_buffsize = 60;
	char* chromo_name;
	char d;
	if ((chromo_name = malloc(chromo_name_buffsize * sizeof(char))) == NULL) {
		printf("Out of memory\n");
		exit(EXIT_FAILURE);
	}
	int chromo_name_index = 0;
	while ((d = getc(ref)) != '\n' && d != ' ') {
		if (chromo_name_index == (chromo_name_buffsize - 1)) {
			chromo_name_buffsize *= 2;
			char* tmp = realloc(chromo_name, chromo_name_buffsize * sizeof(char));
			if (tmp != NULL) {
				chromo_name = tmp;
			}
			else {
				printf("Out of memory\n");
				free(chromo_name);
				exit(EXIT_FAILURE);
			}
		}
		chromo_name[chromo_name_index++] = d;
	}

	chromo_name[chromo_name_index] = '\0';

	return chromo_name;
}

location get_next_loc(FILE* mask_locs) {
	char c;
	char* str;
	if ((str = malloc(100)) == NULL) {
		printf("ERROR: Out of memory\n");
		exit(EXIT_FAILURE);
	}
	size_t buffsize = 100;
	int i = 0;

	while ((c = getc(mask_locs)) != '\n') {
		if (c == EOF) {
			printf("Unexpected end of file\n");
			exit(EXIT_FAILURE);
		}

		/* Need to reallocate memory */
		if (i == (buffsize - 1)) {
			char* tmp = realloc(str, buffsize *= 2);

			if (tmp == NULL) {
				printf("ERROR: Out of memory\n");
				free(str);
				exit(EXIT_FAILURE);
			}

			else {
				str = tmp;
			}
		}

		str[i++] = c;
	}

	str[i] = '\0';
	location loc;

	loc.name = strtok(str, " ");
	loc.start = strtoul(strtok(NULL, " "), NULL, 10);
	loc.end = strtoul(strtok(NULL, " "), NULL, 10);

	free(str);
	
	return loc;

}
