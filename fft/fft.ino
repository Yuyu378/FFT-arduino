// 
// fft.ino
// 
// Feb 13, 2023
// Feb 14, 2023
// 

extern "C" {
  #include <math.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include <stdbool.h>
  
  #include "complex.h"
  #include "numeric.h"
  #include "fft.h"
}



// Functions Declaration --------------------------------------------

// Get remaining memory
int freeRam();
// Print complex type
void cprint(complex a);



// Setup ------------------------------------------------------------
void setup() {
  Serial.begin(9600);
}



// Loop -------------------------------------------------------------
void loop() {

  // Set the length of the data from external sampling
  unsigned N = 32;
  
  // External signal sampling data storage address
  complex* sampled_data = 0;
  if ((sampled_data = (complex*)calloc(N, sizeof(complex))) == NULL) {
		Serial.println(F("Allocation failed!"));
		return 0;
	}

  // Data storage address after initial processing,
  // will be used for the next step of fft function.
  complex* to_transform_data = 0;
  if ((to_transform_data = (complex*)calloc(N, sizeof(complex))) == NULL) {
		Serial.println(F("Allocation failed!"));
		return 0;
	}

  // Data initial processing.
  // After process, store data to the 'to_transform_data' variable
  for (unsigned i = 0; i < N; i++) {
    to_transform_data[i].real = sampled_data[i].real + 5;
    to_transform_data[i].imag = sampled_data[i].imag;
  }
  free(sampled_data);           // has calloc, free it whatever N is

  // Used to store the results of the fft function
  complex* transformed_data = 0;
	transformed_data = fft(to_transform_data, N, backward);
  free(to_transform_data);      // has calloc, free it whatever N is

  // Used to store the data to be inverted
  complex* to_inverse_transform_data = 0;
  if ((to_inverse_transform_data = (complex*)calloc(N, sizeof(complex))) == NULL) {
		Serial.println(F("Allocation failed!"));
		return 0;
	}

  // Data processing.
  // After process, store data to the 'to_inverse_transform_data' variable
  for (unsigned i = 0; i < N; i++) {
    to_inverse_transform_data[i].real = transformed_data[i].real;
    to_inverse_transform_data[i].imag = transformed_data[i].imag;
  }
  if (N > 1) free(transformed_data);             // has not calloc, free it when N > 1

  // Used to store the results of the ifft function
  complex* inverse_transformed_data = ifft(to_inverse_transform_data, N, backward);

  free(to_inverse_transform_data);               // has calloc, free it whatever N is
  if (N > 1) free(inverse_transformed_data);     // has not calloc, free it when N > 1
  
  // Display remaining memory
  Serial.println(freeRam());
  delay(1000);
}



// Functions Implementation -----------------------------------------

int freeRam() {
  extern int __heap_start, *__brkval; 
  int v; 
  return (int) &v - (__brkval == 0 ? (int) &__heap_start : (int) __brkval); 
}

void cprint(complex a) {
	if (a.imag < 0) {
		Serial.print(String(a.real, 6));
		Serial.print(F(" - "));
		Serial.print(String(-a.imag, 6));
		Serial.println(F("i"));
	}
	else {
		Serial.print(String(a.real, 6));
		Serial.print(F(" + "));
		Serial.print(String(a.imag, 6));
		Serial.println(F("i"));
	}
	return;
}
