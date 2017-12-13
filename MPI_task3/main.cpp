#include "matrix.h"
#include "mpi.h"
#include <iostream>
using std::vector;
using std::cout;
using std::endl;

void setTransportVectors(vector<int> &offset, vector<int> &count, const int &state, const int &rowA, const int &rowB, 
	const SparseMatrix &matr1, const SparseMatrix &matr2, vector<int> &linesToProcess)
{
	int ProcRank, ProcNum;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

	switch (state)
	{
	case 0:
	{
		int disp = 0;
		int nextRowPart = rowA;
		int rowCount = matr1.v_row.size();
		for (int i = 0; i < ProcNum; ++i)
		{
			int elCount = 0;
			while (disp < rowCount && matr1.v_row[disp] < nextRowPart)
			{
				++elCount;
				++disp;
			}
			count[i] = elCount;
			if (i != ProcNum - 1)
				offset[i + 1] = offset[i] + elCount;
			nextRowPart += rowA;
		}
	}
		break;
	case 2:
	{
		int disp = 0;
		int nextRowPart = linesToProcess[0];
		int rowCount = matr1.v_row.size();
		for (int i = 0; i < ProcNum; ++i)
		{
			int elCount = 0;
			while (disp < rowCount && matr1.v_row[disp] < nextRowPart)
			{
				++elCount;
				++disp;
			}
			count[i] = elCount;
			if (i != ProcNum - 1)
			{
				offset[i + 1] = offset[i] + elCount;
				nextRowPart += linesToProcess[i + 1];
			}
		}
	}
		break;
	case 1:
	{
		int disp = 0;
		int nextRowPart = rowB;
		int rowCount = matr2.v_row.size();
		for (int i = 0; i < ProcNum; ++i)
		{
			int elCount = 0;
			while (disp < rowCount && matr2.v_row[disp] < nextRowPart)
			{
				++elCount;
				++disp;
			}
			count[i] = elCount;
			if (i != ProcNum - 1)
				offset[i + 1] = offset[i] + elCount;
			nextRowPart += rowB;
		}
	}
		break;
	case 3:
	{
		int disp = 0;
		int nextRowPart = linesToProcess[0];
		int rowCount = matr2.v_row.size();
		for (int i = 0; i < ProcNum; ++i)
		{
			int elCount = 0;
			while (disp < rowCount && matr2.v_row[disp] < nextRowPart)
			{
				++elCount;
				++disp;
			}
			count[i] = elCount;
			if (i != ProcNum - 1)
			{
				offset[i + 1] = offset[i] + elCount;
				nextRowPart += linesToProcess[i + 1];
			}
		}
	}
		break;
	default:
		break;
	}
}


int main(int argc, char ** argv)
{
	int ProcRank, ProcNum;
	double start = 0, end = 0;
	double start2 = 0, end2 = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	int OpToProcess = 0;
	int RowToProcess = 0;
	int ColToProcess = 0;
	int l = 0, m = 0, n = 0;

	//создание типа для передачи структуры Elem
	Elem el;
	MPI_Datatype ELEM_TYPE;
	int block_lengths[2];
	MPI_Aint displacements[2];
	MPI_Aint addresses[3];
	MPI_Datatype typelist[2];
	int blocks_number;
	typelist[0] = MPI_DOUBLE;
	typelist[1] = MPI_DOUBLE;
	block_lengths[0] = block_lengths[1] = 1;
	MPI_Address(&el, &addresses[0]);
	MPI_Address(&(el.Re), &addresses[1]);
	MPI_Address(&(el.Im), &addresses[2]);
	displacements[0] = addresses[1] - addresses[0];
	displacements[1] = addresses[2] - addresses[0];
	blocks_number = 2;
	MPI_Type_struct(blocks_number, block_lengths, displacements, typelist, &ELEM_TYPE);
	MPI_Type_commit(&ELEM_TYPE);

	if (!ProcRank) //введенные размеры матриц
	{
		if (argc > 1)
		{
			l = atoi(argv[1]);
			m = atoi(argv[2]);
			n = atoi(argv[3]);
		}
		else
			std::cin >> l >> m >> n;

		cout << "Matrix A : " << "row = " << l << ", col = " << m << endl;
		cout << "Matrix B : " << "row = " << m << ", col = " << n << endl;
	}

	SparseMatrix matrix(l, m);
	SparseMatrix matrix2(m, n);

	int col = m;

	if (!ProcRank)
	{
		matrix.genMatrix(30);
		matrix2.genMatrix(30);
	}

	SparseMatrix tmatrix(std::move(matrix2.getTransposedV2()));
	SparseMatrix res(l, n);
	SparseMatrix resP(l, n); //результат для параллельной части

	int state = -1;			//номер ситуации
	int rowA = 0;			//количество пересылаемых строк от матрицы matrix
	int rowB = 0;			//количество пересылаемых строк от матрицы tmatrix
	int recElemCountA = 0;	//количество получаемых элементов матрицы 1
	int recElemCountB = 0;	//количество получаемых элементов матрицы 2

	std::vector<int> linesToProcess;

	//transport vectors
	std::vector<int> elemToProcess;
	std::vector<int> offsetToProcess;

	//storage A
	std::vector<Elem> elemStorageA;
	std::vector<int> rowStorageA;
	std::vector<int> colStorageA;
	//storage B
	std::vector<Elem> elemStorageB;
	std::vector<int> rowStorageB;
	std::vector<int> colStorageB;

	//res
	std::vector<Elem> elemPRes;
	std::vector<int> rowPRes;
	std::vector<int> colPRes;

	//results collection
	std::vector<int> resCountFromProc;
	std::vector<int> resOffsetFromProc;
	int sum_of_elems = 0;
	std::vector<Elem> resElemVector;
	std::vector<int> resRowVector;
	std::vector<int> resColVector;

	if (!ProcRank) //последовательная часть
	{		
		start = MPI_Wtime();
		int disp1 = 0; //смещение для матрицы 1
		int matrixSize = matrix.v_row.size();
		for (int row1 = 0; row1 < matrix._row; ++row1)
		{
			int disp2 = 0; //смещение для матрицы 2
			int tmatrixSize = tmatrix.v_row.size();
			for (int row2 = 0; row2 < tmatrix._row; ++row2)
			{
				std::vector<int> col_v(matrix._col, -1);
				Elem sum;
				for (int j = disp1; j < matrixSize && row1 == matrix.v_row[j]; ++j)
				{
					col_v[matrix.v_col[j]] = j;
				}
				for (int j = disp2; j < tmatrixSize && row2 == tmatrix.v_row[j]; ++j)
				{
					if (col_v[tmatrix.v_col[j]] != -1)
						sum += tmatrix.v_elem[j] * matrix.v_elem[col_v[tmatrix.v_col[j]]];
				}
				if (sum != (Elem)0)
				{
					res.v_elem.emplace_back(sum);
					res.v_row.emplace_back(row1);
					res.v_col.emplace_back(row2);
				}
				while (disp2 != tmatrixSize && tmatrix.v_row[disp2] == row2)
					disp2++;
			}
			while (disp1 != matrixSize && matrix.v_row[disp1] == row1)
				disp1++;
		}
		end = MPI_Wtime();

		////матрица 1
		//cout << matrix.v_col.size() << endl; 
		//cout << matrix << endl << endl;
		////матрица 2
		//cout << matrix2.v_col.size() << endl;
		//cout << matrix2 << endl << endl;
		////результат умножения
		//cout << res.v_col.size() << endl;
		//cout << res << endl << endl;

		cout << "Sequential time = " << end - start << endl;
	}	//end

	

	if (!ProcRank) //параллельная часть, используются сгенерированные и транспонированные ранее матрицы (matrix и tmatrix)
	{
		start2 = MPI_Wtime();

		//ProcNum = 2; //отладка

		linesToProcess.resize(ProcNum);

		//transport vectors
		elemToProcess.resize(ProcNum);
		offsetToProcess.resize(ProcNum);

		start2 = MPI_Wtime();
		

		if (l % ProcNum == 0)
		{
			rowA = l / ProcNum;
			rowB = n;
			state = 0;
			setTransportVectors(offsetToProcess, elemToProcess, state, rowA, rowB, matrix, tmatrix, linesToProcess);

			elemStorageB = tmatrix.v_elem;
			rowStorageB = tmatrix.v_row;
			colStorageB = tmatrix.v_col;

			recElemCountB = tmatrix.v_elem.size();
		}
		else if (n % ProcNum == 0)
		{
			rowA = l;
			rowB = n / ProcNum;
			state = 1;
			setTransportVectors(offsetToProcess, elemToProcess, state, rowA, rowB, matrix, tmatrix, linesToProcess);

			elemStorageA = matrix.v_elem;
			rowStorageA = matrix.v_row;
			colStorageA = matrix.v_col;

			recElemCountA = matrix.v_elem.size();
		}
		else if (l > ProcNum)
		{
			int i = 0;
			for (i = 0; i < ProcNum - l % ProcNum; ++i)
				linesToProcess[i] = l / ProcNum;

			for (i; i < ProcNum; ++i)
				linesToProcess[i] = l / ProcNum + 1;
			rowB = n;
			state = 2;
			setTransportVectors(offsetToProcess, elemToProcess, state, rowA, rowB, matrix, tmatrix, linesToProcess);

			elemStorageB = tmatrix.v_elem;
			rowStorageB = tmatrix.v_row;
			colStorageB = tmatrix.v_col;

			recElemCountB = tmatrix.v_elem.size();
		}
		else if (n > ProcNum)
		{
			int i = 0;
			for (i = 0; i < ProcNum - n % ProcNum; ++i)
				linesToProcess[i] = n / ProcNum;

			for (i; i < ProcNum; ++i)
				linesToProcess[i] = n / ProcNum + 1;
			rowA = l;
			state = 3;
			setTransportVectors(offsetToProcess, elemToProcess, state, rowA, rowB, matrix, tmatrix, linesToProcess);

			elemStorageA = matrix.v_elem;
			rowStorageA = matrix.v_row;
			colStorageA = matrix.v_col;

			recElemCountA = matrix.v_elem.size();
		}
		else
		{
			cout << "Bad input value" << endl;
			return -1;
		}


		resCountFromProc.resize(ProcNum);
		resOffsetFromProc.resize(ProcNum);

		/*for (int i = 0; i < tmatrix.v_col.size(); ++i) //отладка
		{
			cout << tmatrix.v_col[i];
		}
		cout << endl;
		for (int i = 0; i < ProcNum; ++i)
		{
			cout << offsetToProcess[i];
		}
		cout << endl;
		for (int i = 0; i < ProcNum; ++i)
		{
			cout << elemToProcess[i];
		}
		cout << endl;
		cout << state;*/
	}


	MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD); //получение количества колонок



	//получение данных
	MPI_Bcast(&state, 1, MPI_INT, 0, MPI_COMM_WORLD);

	switch (state)
	{
	case 0: //+
	{
		MPI_Scatter(elemToProcess.data(), 1, MPI_INT, &recElemCountA, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&recElemCountB, 1, MPI_INT, 0, MPI_COMM_WORLD);

		elemStorageA.resize(recElemCountA);
		rowStorageA.resize(recElemCountA);
		colStorageA.resize(recElemCountA);

		elemStorageB.resize(recElemCountB);
		rowStorageB.resize(recElemCountB);
		colStorageB.resize(recElemCountB);

		MPI_Bcast(elemStorageB.data(), recElemCountB, ELEM_TYPE, 0, MPI_COMM_WORLD);
		MPI_Bcast(rowStorageB.data(), recElemCountB, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(colStorageB.data(), recElemCountB, MPI_INT, 0, MPI_COMM_WORLD);

		MPI_Scatterv(matrix.v_elem.data(), elemToProcess.data(), offsetToProcess.data(), ELEM_TYPE, elemStorageA.data(), recElemCountA, ELEM_TYPE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(matrix.v_row.data(), elemToProcess.data(), offsetToProcess.data(), MPI_INT, rowStorageA.data(), recElemCountA, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(matrix.v_col.data(), elemToProcess.data(), offsetToProcess.data(), MPI_INT, colStorageA.data(), recElemCountA, MPI_INT, 0, MPI_COMM_WORLD);
	}
		break;
	case 2:
	{
		MPI_Scatter(elemToProcess.data(), 1, MPI_INT, &recElemCountA, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&recElemCountB, 1, MPI_INT, 0, MPI_COMM_WORLD);

		elemStorageA.resize(recElemCountA);
		rowStorageA.resize(recElemCountA);
		colStorageA.resize(recElemCountA);

		elemStorageB.resize(recElemCountB);
		rowStorageB.resize(recElemCountB);
		colStorageB.resize(recElemCountB);

		MPI_Bcast(elemStorageB.data(), recElemCountB, ELEM_TYPE, 0, MPI_COMM_WORLD);
		MPI_Bcast(rowStorageB.data(), recElemCountB, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(colStorageB.data(), recElemCountB, MPI_INT, 0, MPI_COMM_WORLD);

		MPI_Scatterv(matrix.v_elem.data(), elemToProcess.data(), offsetToProcess.data(), ELEM_TYPE, elemStorageA.data(), recElemCountA, ELEM_TYPE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(matrix.v_row.data(), elemToProcess.data(), offsetToProcess.data(), MPI_INT, rowStorageA.data(), recElemCountA, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(matrix.v_col.data(), elemToProcess.data(), offsetToProcess.data(), MPI_INT, colStorageA.data(), recElemCountA, MPI_INT, 0, MPI_COMM_WORLD);

	}
		break;
	case 1:
	{
		MPI_Bcast(&recElemCountA, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatter(elemToProcess.data(), 1, MPI_INT, &recElemCountB, 1, MPI_INT, 0, MPI_COMM_WORLD);

		elemStorageA.resize(recElemCountA);
		rowStorageA.resize(recElemCountA);
		colStorageA.resize(recElemCountA);

		elemStorageB.resize(recElemCountB);
		rowStorageB.resize(recElemCountB);
		colStorageB.resize(recElemCountB);

		MPI_Bcast(elemStorageA.data(), recElemCountA, ELEM_TYPE, 0, MPI_COMM_WORLD);
		MPI_Bcast(rowStorageA.data(), recElemCountA, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(colStorageA.data(), recElemCountA, MPI_INT, 0, MPI_COMM_WORLD);

		MPI_Scatterv(tmatrix.v_elem.data(), elemToProcess.data(), offsetToProcess.data(), ELEM_TYPE, elemStorageB.data(), recElemCountB, ELEM_TYPE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(tmatrix.v_row.data(), elemToProcess.data(), offsetToProcess.data(), MPI_INT, rowStorageB.data(), recElemCountB, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(tmatrix.v_col.data(), elemToProcess.data(), offsetToProcess.data(), MPI_INT, colStorageB.data(), recElemCountB, MPI_INT, 0, MPI_COMM_WORLD);
	}
		break;
	case 3:
	{
		MPI_Bcast(&recElemCountA, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatter(elemToProcess.data(), 1, MPI_INT, &recElemCountB, 1, MPI_INT, 0, MPI_COMM_WORLD);

		elemStorageA.resize(recElemCountA);
		rowStorageA.resize(recElemCountA);
		colStorageA.resize(recElemCountA);

		elemStorageB.resize(recElemCountB);
		rowStorageB.resize(recElemCountB);
		colStorageB.resize(recElemCountB);

		MPI_Bcast(elemStorageA.data(), recElemCountA, ELEM_TYPE, 0, MPI_COMM_WORLD);
		MPI_Bcast(rowStorageA.data(), recElemCountA, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(colStorageA.data(), recElemCountA, MPI_INT, 0, MPI_COMM_WORLD);

		MPI_Scatterv(tmatrix.v_elem.data(), elemToProcess.data(), offsetToProcess.data(), ELEM_TYPE, elemStorageB.data(), recElemCountB, ELEM_TYPE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(tmatrix.v_row.data(), elemToProcess.data(), offsetToProcess.data(), MPI_INT, rowStorageB.data(), recElemCountB, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(tmatrix.v_col.data(), elemToProcess.data(), offsetToProcess.data(), MPI_INT, colStorageB.data(), recElemCountB, MPI_INT, 0, MPI_COMM_WORLD);
	}
		break;
	default:
		break;
	}
	//конец получения данных

	//подсчет данных
	int disp1 = 0; //смещение для матрицы 1
	int matrixSize = recElemCountA;
	for (int row1 = rowStorageA[0]; row1 <= rowStorageA[rowStorageA.size() - 1]; ++row1)
	{
		int disp2 = 0; //смещение для матрицы 2
		int tmatrixSize = recElemCountB;
		for (int row2 = rowStorageB[0]; row2 <= rowStorageB[rowStorageB.size() - 1]; ++row2) 
		{
			std::vector<int> col_v(col, -1);
			Elem sum;
			for (int j = disp1; j < matrixSize && row1 == rowStorageA[j]; ++j)
			{
				col_v[colStorageA[j]] = j;
			}
			for (int j = disp2; j < tmatrixSize && row2 == rowStorageB[j]; ++j)
			{
				if (col_v[colStorageB[j]] != -1)
					sum += elemStorageB[j] * elemStorageA[col_v[colStorageB[j]]];
			}
			if (sum != (Elem)0)
			{
				elemPRes.emplace_back(sum);
				rowPRes.emplace_back(row1);
				colPRes.emplace_back(row2);
			}
			while (disp2 != tmatrixSize && rowStorageB[disp2] == row2)
				disp2++;
		}
		while (disp1 != matrixSize && rowStorageA[disp1] == row1)
			disp1++;
	}

	int sendVectorSize = elemPRes.size();
	//подсчет данных

	MPI_Gather(&sendVectorSize, 1, MPI_INT, resCountFromProc.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (!ProcRank)
	{
		for (auto& n : resCountFromProc)
			sum_of_elems += n;
		for (int i = 1; i < ProcNum; ++i)
			resOffsetFromProc[i] = resCountFromProc[i - 1] + resOffsetFromProc[i - 1];

		resElemVector.resize(sum_of_elems);
		resRowVector.resize(sum_of_elems);
		resColVector.resize(sum_of_elems);
	}

	MPI_Gatherv(elemPRes.data(), elemPRes.size(), ELEM_TYPE, resElemVector.data(), resCountFromProc.data(), resOffsetFromProc.data(), ELEM_TYPE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(rowPRes.data(), rowPRes.size(), MPI_INT, resRowVector.data(), resCountFromProc.data(), resOffsetFromProc.data(), MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(colPRes.data(), colPRes.size(), MPI_INT, resColVector.data(), resCountFromProc.data(), resOffsetFromProc.data(), MPI_INT, 0, MPI_COMM_WORLD);

	if (!ProcRank)
	{
		end2 = MPI_Wtime();
		cout << "Parallel time   = " << end2 - start2 << endl;
		cout << "Acceleration    = " << (end - start) / (end2 - start2) << endl;
		resP.v_elem = resElemVector;
		resP.v_row = resRowVector;
		resP.v_col = resColVector;
		cout << "Number of processes = " << ProcNum << endl;
		if (state == 1 || state == 3)
			resP.sortByRow();
		cout << "Matrixes are " << ((resP == res) ? "match" : "different") << endl;
	}

	MPI_Finalize();
	return 0;
}