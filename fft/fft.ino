extern "C" {
  #include <math.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include <stdbool.h>
  
  #include "complex.h"
  #include "numeric.h"
  #include "fft.h"
}

// Set the length of the data from external sampling
unsigned N = 34;

// Data storage address after initial processing,
// will be used for the next step of fft function.
complex* to_transform_data = 0;

// Used to store the results of the fft function
complex* transformed_data = 0;

// Used to store the data to be inverted
complex* to_inverse_transform_data = 0;

// Used to store the results of the ifft function
complex* inverse_transformed_data = 0;

void setup() {
  // put your setup code here, to run once:

  // External signal sampling data storage address
  complex* sampled_data = 0;
  if ((sampled_data = (complex*)calloc(N, sizeof(complex))) == NULL) {
		// Do something if memory allocation failed.
		return 0;
	}
}

// Example
void loop() {
  // put your main code here, to run repeatedly:

  for (unsigned i = 0; i < N; i++){
    // Get the data from the sampling process and store it in the 'sampled_data'
  }

  // Process data in 'sampled_data', like zero-padding/filtering

  // Stroe processed data from 'sampled_data' to 'to_transform_data'

  // Transform data
  transformed_data = fft(to_transform_data, N, backward);

  // Process data in 'transformed_data'

  // Stroe processed data from 'transformed_data' to 'to_inverse_transform_data'

  // Inverse Transform data
  inverse_transformed_data = ifft(to_inverse_transform_data, N, backward);

  // Process data in 'inverse_transformed_data'

  if (N != 1) free(transformed_data);
	if (N != 1) free(inverse_transformed_data);
}
