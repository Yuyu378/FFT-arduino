// 
// fft.ino
// 
//    Main Program
// 
//    v0.1.0
//		----------------------------
// 

extern "C" {
  #include <math.h>
  #include <stdlib.h>
  #include <stdbool.h>
  
  #include "complex.h"
  #include "numeric.h"
  #include "fft.h"
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Functions Declaration
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

// Get remaining memory
int freeRam();
// Print complex type
void cprint(_Dcomplex a);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Setup
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void setup() {
  Serial.begin(9600);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Main Loop
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void loop() {

  unsigned long StartingTime, EndingTime, ElapsedMicroseconds;

  // Set the length of the data from external sampling
  unsigned N = 119;
  bool isprint = true;
  
  _Dcomplex* x = 0;
	if ((x = (_Dcomplex*)calloc(N, sizeof(_Dcomplex))) == NULL) {
		Serial.println(F("Allocation failed!"));
		return 0;
	}
  
  Serial.println(F("Raw data :"));
  StartingTime = micros();

  for (unsigned i = 0; i < N; i++) {
		x[i] = _Cbuild(i, 0.);
		if (isprint) cprint(x[i]);
	}

  EndingTime = micros();
  ElapsedMicroseconds = EndingTime - StartingTime;
  Serial.println("run_time: " + String(ElapsedMicroseconds) + " \u03bcs\n");

  {
    Serial.println(F("Transformed data :"));
    StartingTime = micros();
    
    _Dcomplex* y = fft(x, N, backward);
    
    EndingTime = micros();
    ElapsedMicroseconds = EndingTime - StartingTime;
    if (isprint) {
			for (unsigned i = 0; i < N; i++) {
				cprint(*(y + i));
			}
		}
		Serial.println("run_time: " + String(ElapsedMicroseconds) + " \u03bcs\n");
    if (N != 1) free(y);
  }
  
  free(x);
  
  // Display remaining memory
  Serial.println(freeRam());
  Serial.println();
  delay(1000);
  //while(1) {}            // Stop loop()
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Functions Implementation
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

int freeRam() {
  extern int __heap_start, *__brkval; 
  int v; 
  return (int) &v - (__brkval == 0 ? (int) &__heap_start : (int) __brkval); 
}

void cprint(_Dcomplex a) {
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
