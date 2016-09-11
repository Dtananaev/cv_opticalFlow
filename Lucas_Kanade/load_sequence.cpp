#include "load_sequence.h"

CVector< CMatrix<float> > loadSequence( std::string Filename)
{

	CVector<std::string> Input;
	std::string InputDir;

	// Determine input directory
	std::string s = Filename;
	s.erase(s.find_last_of("/") + 1, s.length());
	InputDir = s;
	s = Filename;
	s.erase(0, s.find_last_of('.'));

	if (s == ".txt" || s == ".TXT") {
		int aImageCount;
        std::ifstream aStream(Filename.c_str());
		aStream >> aImageCount;
		Input.setSize(aImageCount);
		for (int i = 0; i < aImageCount; i++) {
			std::string s;
			aStream >> s;
			Input(i) = InputDir + s;
		}
	}
	else {
		std::cout << "Must pass a txt file as input" << std::endl;
		exit(1);
	}

    CVector< CMatrix<float> > noisy_images;
	noisy_images.setSize(Input.size());

	for (int i = 0; i < Input.size(); i++)
	{
		noisy_images(i).readFromPGM(Input(i).c_str());
	}
	return noisy_images;
}

