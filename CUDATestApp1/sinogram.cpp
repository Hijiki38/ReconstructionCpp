#include"sinogram.h"
#include"methods.h"


namespace Reconstruction {

	VectorXf sinogram::get_sinovec() { return eigen_sinovec; }
	
	int sinogram::get_nd() { return n_d; }
	
	int sinogram::get_nv() { return n_v; }

	sinogram* sinogram::read_sinogram(string inpath) {

		//ifstream stream(inpath);
		//FILE* fp;

		//if ((fopen_s(&fp, inpath.c_str(), "r")) != 0) {
		//	std::cout << "File not found.";
		//	return nullptr;
		//}

		//size_t buffersize = 5 * 1024 * 1024;
		//char* buffer_infile = static_cast<char*>(malloc(buffersize));
		//char buffer_value[1024];

		//string line;
		//vector<float> *sinovec;
		//int count_row = 0;
		//int count_all = 0;

		//cout << "read sinogram";

		//int nrow = 0;
		//int nelem = 0;

		//size_t size;

		//while (!feof(fp)) 
		//{
		//	size = fread(&buffer_infile[0], 1, buffersize / sizeof(char), fp);
		//	for (int i = 0; i < size; i++)
		//	{
		//		if (buffer_infile[i] == ',')
		//		{
		//			nelem++;
		//		}
		//		else if (buffer_infile[i] == '\n')
		//		{
		//			if (buffer_infile[i - 1] != ',')
		//			{
		//				nelem++;
		//			}
		//			nrow++;
		//		}
		//	}
		//}

		//sinovec = new vector<float>(nelem);

		//cout << "\nsizebuf:" << buffersize;
		//cout << "\nsize:" << size;

		//for (int i = 0, counter = 0; i < size; i++) {
		//	if ((buffer_infile[i] >= (char)'0' && buffer_infile[i] <= (char)'9') 
		//		|| buffer_infile[i] == (char)'.' || buffer_infile[i] == (char)'-')
		//	{
		//		buffer_value[counter++] = buffer_infile[i];
		//	}
		//	else if (buffer_infile[i] == (char)'\n' || buffer_infile[i] == (char)',') 
		//	{
		//		buffer_value[counter] = '\0';
		//		sinovec->push_back(atof(buffer_value));
		//		cout << "\nsinovec:" << atof(buffer_value);
		//		counter = 0;
		//	}
		//}

		//free(buffer_infile);

		//cout << "\n uoooooo";
		//cout << "\nnrow:" << nrow << ", nelem:" << nelem;
		//cout << "\nsize of sinogram: (" << nrow << ", " << nelem / nrow << ")";

		//sinogram* sino = new sinogram(sinovec, nelem / nrow, nrow);

		//cout << "\nd: " << (*sino).get_nd() << ", v: " << (*sino).get_nv();

		//return sino;

		//while (getline(stream, line)) {
		//	cout << "\rreadrow:" << count_row;
		//	vector<string> strs = Reconstruction::splitstring(line, ',');
		//	for (int i = 0; i < strs.size(); i++) {
		//		sinovec->push_back(stof(strs.at(i)));
		//		count_all++;
		//		cout << "\rreadelem:" << count_all;
		//	}
		//	count_row++;
		//}

		//cout << "\nd: " << (count_all / count_row) << ", v: " << count_row;

		//sinogram* sino = new sinogram(sinovec, count_all / count_row, count_row);

		//cout << "\nd: " << (*sino).get_nd() << ", v: " << (*sino).get_nv();

		//return sino;


		ifstream stream(inpath);
		string line;
		vector<float> sinovec;
		int count_row = 0;
		int count_all = 0;

		cout << "read sinogram";

		while (getline(stream, line)) {
			cout << "readrow:" << count_row << "\n";
			vector<string> strs = Reconstruction::splitstring(line, ',');
			for (int i = 0; i < strs.size(); i++) {
				sinovec.push_back(stof(strs.at(i)));
				count_all++;
				cout << "readelem:" << count_all << "\n";
			}
			count_row++;
		}

		cout << "d: " << (count_all / count_row) << ", v: " << count_row << "\n";

		sinogram* sino = new sinogram(&sinovec, count_all / count_row, count_row);

		cout << "d: " << (*sino).get_nd() << ", v: " << (*sino).get_nv() << "\n";

		return sino;
	}
}