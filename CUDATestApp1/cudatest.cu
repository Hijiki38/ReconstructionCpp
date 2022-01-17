#include <stdio.h>
#include <iostream>
#include <time.h>

// GPU�Ōv�Z����ۂ̊֐�
__global__ void gpu_function(float* d_x, float* d_y)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	d_y[i] = sin(d_x[i]) * sin(d_x[i]) + cos(d_x[i]) * cos(d_x[i]);
}

// CPU�Ōv�Z����ۂ̊֐�
void cpu_function(int n, float* x, float* y)
{
	for (int i = 0; i < n; i++) {
		y[i] = sin(x[i]) * sin(x[i]) + cos(x[i]) * cos(x[i]);
	}
}

// main function
int cudamain(void)
{
	bool GPU = true;

	int N = 1000000;
	float* host_x, * host_y, * dev_x, * dev_y;

	// CPU���̗̈�m��
	host_x = (float*)malloc(N * sizeof(float));
	host_y = (float*)malloc(N * sizeof(float));

	// �����l����͂���
	for (int i = 0; i < N; i++) {
		host_x[i] = rand();
	}

	int start = clock();

	if (GPU == true) {

		// �f�o�C�X(GPU)���̗̈�m��
		cudaMalloc(&dev_x, N * sizeof(float));
		cudaMalloc(&dev_y, N * sizeof(float));

		// CPU��GPU�̃f�[�^�R�s�[
		cudaMemcpy(dev_x, host_x, N * sizeof(float), cudaMemcpyHostToDevice);

		// GPU�Ōv�Z
		gpu_function << <(N + 255) / 256, 256 >> > (dev_x, dev_y);

		// GPU��CPU�̃f�[�^�R�s�[
		cudaMemcpy(host_y, dev_y, N * sizeof(float), cudaMemcpyDeviceToHost);

	}
	else {
		// CPU�Ōv�Z
		cpu_function(N, host_x, host_y);
	}

	int end = clock();

	// �v�Z���������s���Ă��邩�m�F
	float sum = 0.0f;
	for (int j = 0; j < N; j++) {
		sum += host_y[j];
	}
	std::cout << sum << std::endl;

	// �Ō�Ɍv�Z���Ԃ�\��
	std::cout << end - start << "[ms]" << std::endl;

	return 0;
}