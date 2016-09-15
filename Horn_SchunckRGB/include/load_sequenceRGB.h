#ifndef LOAD_SEQUENCERGB_H_
#define LOAD_SEQUENCERGB_H_
#include <CVector.h>
#include <CMatrix.h>
#include <CTensor.h>
#include <string>

CVector< CTensor<float> > loadSequence(std::string Filename);

#endif
