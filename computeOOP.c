// Read dump file and data file
// Read the config file containing information about vectors
// Check the data file for bond connection.
// Draw a vector only if those two atoms are covalently connected.
// Save info about all the vectors
// Then calculate orientation order parameter

// This code should be written in a generic manner, such that it'll work for all systems

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>

typedef struct orderParameter
{
	int atom1, atom2, atom3, atom4;
	float distance, theta_rad, theta_deg, orderParameter;
} ORDERPARAMETER;

typedef struct dumpinfo
{
	int timestep, nAtoms;
	float xlo, xhi, ylo, yhi, zlo, zhi;
} DUMPFILE_INFO;

typedef struct config
{
	int atom1, atom2;
} CONFIG;

typedef struct datafile_atoms
{
	int resNumber;
	char resName[6], atomName[6], atomType2[6], molName[6];

	int id, molType, atomType;
	float charge, x, y, z;
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2, atom1Type, atom2Type;
} DATA_BONDS;

typedef struct datafile_angles
{
	int id, angleType, atom1, atom2, atom3;
} DATA_ANGLES;

typedef struct datafile_dihedrals
{
	int id, dihedralType, atom1, atom2, atom3, atom4;
} DATA_DIHEDRALS;

typedef struct datafile_impropers
{
	int id, improperType, atom1, atom2, atom3, atom4;
} DATA_IMPROPERS;

typedef struct datafileInfo
{
	int nAtoms, nBonds, nAngles, nDihedrals, nImpropers;
	int nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes;
} DATAFILE_INFO;

char *getInputFileName()
{
	int nFiles = 0;
	char *inputFileName, fileExtension[200]/*, terminalString[200]*/;
	int fileRequired;
	inputFileName = (char *) malloc (200 * sizeof (char));

	getFilenameAgain:
	printf("Enter the file extension or a match string to search in current directory...\n --> "); scanf ("%s", fileExtension);
	// fgets (terminalString, sizeof (terminalString), stdin);
	// sscanf (terminalString, "%s", fileExtension); 
	printf("\n");
	nFiles = displayFiles(fileExtension);

	if (nFiles > 0)
	{
		fprintf(stdout, "\nWhich file would you like to input? Enter a number between (1 to %d): ", nFiles); 
		fflush (stdout);
		scanf ("%d", &fileRequired);
		// fgets(terminalString, sizeof (terminalString), stdin);
		// sscanf (terminalString, "%d", &fileRequired); 
	}
	else
	{
		printf("No files found with the match string. Try again!\n"); goto getFilenameAgain;
	}

	nFiles = 0;
	DIR *parentDirectory;
	parentDirectory = opendir ("./");

	struct dirent *filePointer;

	/* Scan all the files using filePointer */
	while ((filePointer = readdir (parentDirectory)))
	{
		if (isFile(filePointer -> d_name) && strstr(filePointer -> d_name, fileExtension))
		{
			nFiles++;
			if (fileRequired == nFiles)
			{
				strcpy (inputFileName, filePointer -> d_name);
			}
		}
	}
	return inputFileName;
}

DATAFILE_INFO readData (FILE *input, DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, DATA_IMPROPERS **impropers)
{
	printf("Reading LAMMPS data file...\n");
	rewind (input);

	int isAtomLine = 0, /*nAtoms = -1,*/ nAtomLine = 0;
	int isBondLine = 0, /*nBonds = -1,*/ nBondLine = 0;
	int isAngleLine = 0, /*nAngles = -1,*/ nAngleLine = 0;
	int isDihedralLine = 0, /*nDihedrals = -1,*/ nDihedralLine = 0;
	int isImproperLine = 0, /*nImpropers = -1,*/ nImproperLine = 0;
	int printHeaderInfo = 1;

	DATAFILE_INFO datafile;
	datafile.nAtoms = -1;
	datafile.nBonds = -1;
	datafile.nAngles = -1;
	datafile.nDihedrals = -1;
	datafile.nImpropers = -1;

	char lineString[1000];

	*atoms = NULL;
	*bonds = NULL;
	*angles = NULL;
	*dihedrals = NULL;
	*impropers = NULL;

	while ((fgets (lineString, 1000, input) != NULL))
	{
		if (strstr (lineString, "atoms"))
		{
			sscanf (lineString, "%d \n", &datafile.nAtoms);
			fprintf(stdout, "nAtoms detected: %d\n", datafile.nAtoms);
			fflush (stdout);
			(*atoms) = (DATA_ATOMS *) malloc (datafile.nAtoms * sizeof (DATA_ATOMS));
		}

		if (strstr (lineString, "bonds"))
		{
			sscanf (lineString, "%d \n", &datafile.nBonds);
			fprintf(stdout, "nBonds detected: %d\n", datafile.nBonds);
			fflush (stdout);
			(*bonds) = (DATA_BONDS *) malloc (datafile.nBonds * sizeof (DATA_BONDS));
		}

		if (strstr (lineString, "angles"))
		{
			sscanf (lineString, "%d \n", &datafile.nAngles);
			fprintf(stdout, "nAngles detected: %d\n", datafile.nAngles);
			fflush (stdout);
			(*angles) = (DATA_ANGLES *) malloc (datafile.nAngles * sizeof (DATA_ANGLES));
		}

		if (strstr (lineString, "dihedrals"))
		{
			sscanf (lineString, "%d \n", &datafile.nDihedrals);
			fprintf(stdout, "nDihedrals detected: %d\n", datafile.nDihedrals);
			fflush (stdout);
			(*dihedrals) = (DATA_DIHEDRALS *) malloc (datafile.nDihedrals * sizeof (DATA_DIHEDRALS));
		}

		if (strstr (lineString, "impropers"))
		{
			sscanf (lineString, "%d \n", &datafile.nImpropers);
			fprintf(stdout, "nImpropers detected: %d\n", datafile.nImpropers);
			fflush (stdout);
			(*impropers) = (DATA_IMPROPERS *) malloc (datafile.nImpropers * sizeof (DATA_IMPROPERS));
		}

		if (strstr (lineString, "atom types"))
		{
			sscanf (lineString, "%d \n", &datafile.nAtomTypes);
			fprintf(stdout, "%s%d\n", "nAtomTypes detected: ", datafile.nAtomTypes);
			fflush (stdout);
		}

		if (strstr (lineString, "bond types"))
		{
			sscanf (lineString, "%d \n", &datafile.nBondTypes);
			fprintf(stdout, "%s%d\n", "nBondTypes detected: ", datafile.nBondTypes);
			fflush (stdout);
		}

		if (strstr (lineString, "angle types"))
		{
			sscanf (lineString, "%d \n", &datafile.nAngleTypes);
			fprintf(stdout, "%s%d\n", "nAngleTypes detected: ", datafile.nAngleTypes);
			fflush (stdout);
		}

		if (strstr (lineString, "dihedral types"))
		{
			sscanf (lineString, "%d \n", &datafile.nDihedralTypes);
			fprintf(stdout, "%s%d\n", "nDihedralTypes detected: ", datafile.nDihedralTypes);
			fflush (stdout);
		}

		if (strstr (lineString, "improper types"))
		{
			sscanf (lineString, "%d \n", &datafile.nImproperTypes);
			fprintf(stdout, "%s%d\n", "nImproperTypes detected: ", datafile.nImproperTypes);
			fflush (stdout);
		}

		if ((datafile.nAtoms >= 0) && (datafile.nBonds >= 0) && (datafile.nAngles >= 0) && (datafile.nDihedrals >= 0) && (datafile.nImpropers >= 0) && (printHeaderInfo))
			printHeaderInfo = 0;

		if (strstr (lineString, "Atoms"))
		{
			isAtomLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Bonds"))
		{
			isBondLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Angles"))
		{
			isAngleLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Dihedrals"))
		{
			isDihedralLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Impropers"))
		{
			isImproperLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (isAtomLine)
		{
			sscanf (lineString, "%d %d %d %f %f %f %f\n", 
				&(*atoms)[nAtomLine].id, 
				&(*atoms)[nAtomLine].molType, 
				&(*atoms)[nAtomLine].atomType, 
				&(*atoms)[nAtomLine].charge, 
				&(*atoms)[nAtomLine].x, 
				&(*atoms)[nAtomLine].y, 
				&(*atoms)[nAtomLine].z);
			nAtomLine++;
			if (nAtomLine == datafile.nAtoms)
				isAtomLine = 0;
		}

		if (isBondLine)
		{
			sscanf (lineString, "%d %d %d %d\n", 
				&(*bonds)[nBondLine].id, 
				&(*bonds)[nBondLine].bondType, 
				&(*bonds)[nBondLine].atom1, 
				&(*bonds)[nBondLine].atom2);
			nBondLine++;
			if (nBondLine == datafile.nBonds)
				isBondLine = 0;
		}

		if (isAngleLine)
		{
			sscanf (lineString, "%d %d %d %d %d\n", 
				&(*angles)[nAngleLine].id, 
				&(*angles)[nAngleLine].angleType, 
				&(*angles)[nAngleLine].atom1, 
				&(*angles)[nAngleLine].atom2, 
				&(*angles)[nAngleLine].atom3);
			nAngleLine++;
			if (nAngleLine == datafile.nAngles)
				isAngleLine = 0;
		}

		if (isDihedralLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*dihedrals)[nDihedralLine].id, 
				&(*dihedrals)[nDihedralLine].dihedralType, 
				&(*dihedrals)[nDihedralLine].atom1, 
				&(*dihedrals)[nDihedralLine].atom2, 
				&(*dihedrals)[nDihedralLine].atom3, 
				&(*dihedrals)[nDihedralLine].atom4);
			nDihedralLine++;
			if (nDihedralLine == datafile.nDihedrals)
				isDihedralLine = 0;
		}

		if (isImproperLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*impropers)[nImproperLine].id, 
				&(*impropers)[nImproperLine].improperType, 
				&(*impropers)[nImproperLine].atom1, 
				&(*impropers)[nImproperLine].atom2, 
				&(*impropers)[nImproperLine].atom3, 
				&(*impropers)[nImproperLine].atom4);
			nImproperLine++;
			if (nImproperLine == datafile.nImpropers)
				isImproperLine = 0;
		}
	}

	printf("\nFrom input data file:\n\n nAtoms: %d\n nBonds: %d\n nAngles: %d\n nDihedrals: %d\n nImpropers: %d\n\n", datafile.nAtoms, datafile.nBonds, datafile.nAngles, datafile.nDihedrals, datafile.nImpropers);

	rewind (input);
	return datafile;
}

int isFile(const char *name)
{
	DIR *directory = opendir (name);
	if (directory!=NULL)
	{
		closedir(directory);
		return 0;
	}
	if(errno==ENOTDIR)
	{
		return 1;
	}

	return -1;
}

int displayFiles(const char *fileExtension)
{
	int nFiles = 0;
	DIR *parentDirectory;
	parentDirectory = opendir ("./");

	struct dirent *filePointer;
	/* Scan all the files using filePointer */
	while ((filePointer = readdir (parentDirectory)))
	{
		if (isFile(filePointer -> d_name) && strstr(filePointer -> d_name, fileExtension))
		{
			nFiles++;
			printf("%d --> %s\n", nFiles, filePointer -> d_name);
		}
	}
	return nFiles;
}

CONFIG *readConfig (FILE *inputConfigFile)
{
	rewind (inputConfigFile);
	CONFIG *inputVectors;
	char lineString[1000];
	int nLines = 0;

	while (fgets (lineString, 1000, inputConfigFile) != NULL)
	{
		if (lineString[0] != '#')
		{
			nLines++;
		}
	}

	inputVectors = (CONFIG *) malloc (nLines * sizeof (CONFIG));
	rewind (inputConfigFile);
	nLines = 0;

	while (fgets (lineString, 1000, inputConfigFile) != NULL)
	{
		if (lineString[0] != '#')
		{
			sscanf (lineString, "%d %d\n", &inputVectors[nLines].atom1, &inputVectors[nLines].atom2);
			nLines++;
		}
	}

	rewind (inputConfigFile);
	return inputVectors;
}

DUMPFILE_INFO getDumpFileInfo (FILE *inputDumpFile)
{
	rewind (inputDumpFile);
	char lineString[1000];
	DUMPFILE_INFO dumpfile;

	for (int i = 0; i < 9; ++i)
	{
		fgets (lineString, 1000, inputDumpFile);
		if (i == 1)
		{
			sscanf (lineString, "%d", &dumpfile.timestep);
		}
		if (i == 3)
		{
			sscanf (lineString, "%d", &dumpfile.nAtoms);
		}
		if (i == 5)
		{
			sscanf (lineString, "%f %f\n", &dumpfile.xlo, &dumpfile.xhi);
		}
		if (i == 6)
		{
			sscanf (lineString, "%f %f\n", &dumpfile.ylo, &dumpfile.yhi);
		}
		if (i == 7)
		{
			sscanf (lineString, "%f %f\n", &dumpfile.zlo, &dumpfile.zhi);
		}
	}

	rewind (inputDumpFile);
	return dumpfile;
}

ORDERPARAMETER *printOrderParameter (DATA_ATOMS *dumpAtoms, DUMPFILE_INFO dumpfile, DATAFILE_INFO datafile, DATA_BONDS *bonds, CONFIG *inputVectors, int currentTimestep, unsigned int nElements)
{
	FILE *allData;
	char *allData_string;
	allData_string = (char *) malloc (50 * sizeof (char));
	sprintf (allData_string, "logs/allData_%d.oop", currentTimestep);
	allData = fopen (allData_string, "w");

	ORDERPARAMETER *allData_array;
	allData_array = (ORDERPARAMETER *) malloc (nElements * sizeof (ORDERPARAMETER));

	unsigned int currentElement = 0;

	float x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, distance, dotProduct, magnitude1, magnitude2, cosTheta, theta, orderParameter;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		bonds[i].atom1Type = dumpAtoms[bonds[i].atom1 - 1].atomType;
		bonds[i].atom2Type = dumpAtoms[bonds[i].atom2 - 1].atomType;

		// Checking if the bonds correspond to inputVectors[0]; from the first line of the config file
		if ((bonds[i].atom1Type == inputVectors[0].atom1 && bonds[i].atom2Type == inputVectors[0].atom2) || (bonds[i].atom1Type == inputVectors[0].atom2 && bonds[i].atom2Type == inputVectors[0].atom1))
		{
			for (int j = (i + 1); j < datafile.nBonds; ++j)
			{
				bonds[j].atom1Type = dumpAtoms[bonds[j].atom1 - 1].atomType;
				bonds[j].atom2Type = dumpAtoms[bonds[j].atom2 - 1].atomType;

				if ((bonds[j].atom1Type == inputVectors[1].atom1 && bonds[j].atom2Type == inputVectors[1].atom2) || (bonds[j].atom1Type == inputVectors[1].atom2 && bonds[j].atom2Type == inputVectors[1].atom1))
				{
					// Finding the center of two bonds (x1, y1, z1) and (x2, y2, z2)
					x1 = (dumpAtoms[bonds[i].atom1 - 1].x + dumpAtoms[bonds[i].atom2 - 1].x) / 2; x2 = (dumpAtoms[bonds[j].atom1 - 1].x + dumpAtoms[bonds[j].atom2 - 1].x) / 2; y1 = (dumpAtoms[bonds[i].atom1 - 1].y + dumpAtoms[bonds[i].atom2 - 1].y) / 2; y2 = (dumpAtoms[bonds[j].atom1 - 1].y + dumpAtoms[bonds[j].atom2 - 1].y) / 2; z1 = (dumpAtoms[bonds[i].atom1 - 1].z + dumpAtoms[bonds[i].atom2 - 1].z) / 2; z2 = (dumpAtoms[bonds[j].atom1 - 1].z + dumpAtoms[bonds[j].atom2 - 1].z) / 2;

					// Distance between the centers of two bonds
					distance = sqrt (((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) + ((z2 - z1) * (z2  - z1)));

					// Storing the positions of all 4 atoms forming the two bonds of interest
					x1 = dumpAtoms[bonds[i].atom1 - 1].x; y1 = dumpAtoms[bonds[i].atom1 - 1].y; z2 = dumpAtoms[bonds[i].atom1 - 1].z; x2 = dumpAtoms[bonds[i].atom2 - 1].x; y2 = dumpAtoms[bonds[i].atom2 - 1].y; z2 = dumpAtoms[bonds[i].atom2 - 1].z; x3 = dumpAtoms[bonds[j].atom1 - 1].x; y3 = dumpAtoms[bonds[j].atom1 - 1].y; z3 = dumpAtoms[bonds[j].atom1 - 1].z; x4 = dumpAtoms[bonds[j].atom2 - 1].x; y4 = dumpAtoms[bonds[j].atom2 - 1].y; z4 = dumpAtoms[bonds[j].atom2 - 1].z; dotProduct = ((x2 - x1) * (x4 - x3)) + ((y2 - y1) * (y4 - y3)) + ((z2 - z1) * (z4 - z3)); magnitude1 = ((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) + ((z2 - z1) * (z2 - z1)); magnitude2 = ((x4 - x3) * (x4 - x3)) + ((y4 - y3) * (y4 - y3)) + ((z4 - z3) * (z4 - z3)); cosTheta = dotProduct / (sqrt (magnitude1) * sqrt (magnitude2)); theta = acosf (cosTheta); orderParameter = ((3 * cosTheta * cosTheta) - 1) / 2;

					fprintf(allData, "%d %d %d %d %f %f %f %f\n", bonds[i].atom1, bonds[i].atom2, bonds[j].atom1, bonds[j].atom2, distance, theta, theta * 57.2958, orderParameter);

					allData_array[currentElement].atom1 = bonds[i].atom1; allData_array[currentElement].atom2 = bonds[i].atom2; allData_array[currentElement].atom3 = bonds[j].atom1; allData_array[currentElement].atom4 = bonds[j].atom2; allData_array[currentElement].distance = distance; allData_array[currentElement].theta_rad = theta; allData_array[currentElement].theta_deg = theta * 57.2958; allData_array[currentElement].orderParameter = orderParameter; currentElement++;
				}
			}
		}

		// Checking if the bonds correspond to inputVectors[1]; from the second line of the config file
		else if ((bonds[i].atom1Type == inputVectors[1].atom1 && bonds[i].atom2Type == inputVectors[1].atom2) || (bonds[i].atom1Type == inputVectors[1].atom2 && bonds[i].atom2Type == inputVectors[1].atom1))
		{
			for (int j = (i + 1); j < datafile.nBonds; ++j)
			{
				bonds[j].atom1Type = dumpAtoms[bonds[j].atom1 - 1].atomType;
				bonds[j].atom2Type = dumpAtoms[bonds[j].atom2 - 1].atomType;

				if ((bonds[j].atom1Type == inputVectors[0].atom1 && bonds[j].atom2Type == inputVectors[0].atom2) || (bonds[j].atom1Type == inputVectors[0].atom2 && bonds[j].atom2Type == inputVectors[0].atom1))
				{
					// Finding the center of two bonds (x1, y1, z1) and (x2, y2, z2)
					x1 = (dumpAtoms[bonds[i].atom1 - 1].x + dumpAtoms[bonds[i].atom2 - 1].x) / 2; x2 = (dumpAtoms[bonds[j].atom1 - 1].x + dumpAtoms[bonds[j].atom2 - 1].x) / 2; y1 = (dumpAtoms[bonds[i].atom1 - 1].y + dumpAtoms[bonds[i].atom2 - 1].y) / 2; y2 = (dumpAtoms[bonds[j].atom1 - 1].y + dumpAtoms[bonds[j].atom2 - 1].y) / 2; z1 = (dumpAtoms[bonds[i].atom1 - 1].z + dumpAtoms[bonds[i].atom2 - 1].z) / 2; z2 = (dumpAtoms[bonds[j].atom1 - 1].z + dumpAtoms[bonds[j].atom2 - 1].z) / 2;

					// Distance between the centers of two bonds
					distance = sqrt (((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) + ((z2 - z1) * (z2  - z1)));

					// Storing the positions of all 4 atoms forming the two bonds of interest
					x1 = dumpAtoms[bonds[i].atom1 - 1].x; y1 = dumpAtoms[bonds[i].atom1 - 1].y; z2 = dumpAtoms[bonds[i].atom1 - 1].z; x2 = dumpAtoms[bonds[i].atom2 - 1].x; y2 = dumpAtoms[bonds[i].atom2 - 1].y; z2 = dumpAtoms[bonds[i].atom2 - 1].z; x3 = dumpAtoms[bonds[j].atom1 - 1].x; y3 = dumpAtoms[bonds[j].atom1 - 1].y; z3 = dumpAtoms[bonds[j].atom1 - 1].z; x4 = dumpAtoms[bonds[j].atom2 - 1].x; y4 = dumpAtoms[bonds[j].atom2 - 1].y; z4 = dumpAtoms[bonds[j].atom2 - 1].z; dotProduct = ((x2 - x1) * (x4 - x3)) + ((y2 - y1) * (y4 - y3)) + ((z2 - z1) * (z4 - z3)); magnitude1 = ((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) + ((z2 - z1) * (z2 - z1)); magnitude2 = ((x4 - x3) * (x4 - x3)) + ((y4 - y3) * (y4 - y3)) + ((z4 - z3) * (z4 - z3)); cosTheta = dotProduct / (sqrt (magnitude1) * sqrt (magnitude2)); theta = acosf (cosTheta); orderParameter = ((3 * cosTheta * cosTheta) - 1) / 2;

					fprintf(allData, "%d %d %d %d %f %f %f %f\n", bonds[i].atom1, bonds[i].atom2, bonds[j].atom1, bonds[j].atom2, distance, theta, theta * 57.2958, orderParameter);

					allData_array[currentElement].atom1 = bonds[i].atom1; allData_array[currentElement].atom2 = bonds[i].atom2; allData_array[currentElement].atom3 = bonds[j].atom1; allData_array[currentElement].atom4 = bonds[j].atom2; allData_array[currentElement].distance = distance; allData_array[currentElement].theta_rad = theta; allData_array[currentElement].theta_deg = theta * 57.2958; allData_array[currentElement].orderParameter = orderParameter; currentElement++;
				}				
			}
		}
	}

	fclose (allData);
	// free (allData_array);
	return allData_array;
}

unsigned int getNElements (DATAFILE_INFO datafile, DATA_ATOMS *dumpAtoms, DATA_BONDS *bonds, CONFIG *inputVectors)
{
	unsigned int nElements = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		bonds[i].atom1Type = (int) dumpAtoms[bonds[i].atom1 - 1].atomType;
		bonds[i].atom2Type = (int) dumpAtoms[bonds[i].atom2 - 1].atomType;

		if ((bonds[i].atom1Type == inputVectors[0].atom1 && bonds[i].atom2Type == inputVectors[0].atom2) || (bonds[i].atom1Type == inputVectors[0].atom2 && bonds[i].atom2Type == inputVectors[0].atom1))
		{
			for (int j = (i + 1); j < datafile.nBonds; ++j)
			{
				bonds[j].atom1Type = dumpAtoms[bonds[j].atom1 - 1].atomType;
				bonds[j].atom2Type = dumpAtoms[bonds[j].atom2 - 1].atomType;

				if ((bonds[j].atom1Type == inputVectors[1].atom1 && bonds[j].atom2Type == inputVectors[1].atom2) || (bonds[j].atom1Type == inputVectors[1].atom2 && bonds[j].atom2Type == inputVectors[1].atom1))
				{
					nElements++;
				}
			}
		}
		else if ((bonds[i].atom1Type == inputVectors[1].atom1 && bonds[i].atom2Type == inputVectors[1].atom2) || (bonds[i].atom1Type == inputVectors[1].atom2 && bonds[i].atom2Type == inputVectors[1].atom1))
		{
			for (int j = (i + 1); j < datafile.nBonds; ++j)
			{
				bonds[j].atom1Type = dumpAtoms[bonds[j].atom1 - 1].atomType;
				bonds[j].atom2Type = dumpAtoms[bonds[j].atom2 - 1].atomType;

				if ((bonds[j].atom1Type == inputVectors[0].atom1 && bonds[j].atom2Type == inputVectors[0].atom2) || (bonds[j].atom1Type == inputVectors[0].atom2 && bonds[j].atom2Type == inputVectors[0].atom1))
				{
					nElements++;
				}				
			}
		}

	}

	return nElements;
}

void computeOrderParameter (FILE *inputDumpFile, DATAFILE_INFO datafile, DATA_BONDS *bonds, CONFIG *inputVectors)
{
	DUMPFILE_INFO dumpfile;
	dumpfile = getDumpFileInfo (inputDumpFile);

	char lineString[1000];
	int isTimestep = 0, currentTimestep, currentLine = 0, isFirstTimestep = 1, currentDumpstep = 0;

	unsigned int nElements = 0;
	int isNElementsSet = 0;

	// Datafile struct is used to store dump atom information
	DATA_ATOMS *dumpAtoms;
	dumpAtoms = (DATA_ATOMS *) malloc (dumpfile.nAtoms * sizeof (DATA_ATOMS));

	ORDERPARAMETER *allData_array;

	printf("\n");

	// Reading and processing dump information
	while (fgets (lineString, 1000, inputDumpFile) != NULL)
	{
		if (currentLine == (9 + dumpfile.nAtoms) && nElements == 0)
		{
			nElements = getNElements (datafile, dumpAtoms, bonds, inputVectors);
			printf("Allocating for %lu elements...\n", nElements);
			allData_array = (ORDERPARAMETER *) malloc (nElements * sizeof (ORDERPARAMETER));
			printf("Memory allocated successfully...\n");
		}

		if (currentDumpstep > 2 && nElements > 0)
		{
			sscanf (lineString, "%d", &currentTimestep);
			printf("Scanning timestep: %d...\n", currentTimestep);
			fflush (stdout); 

			printf("nElements: %lu\n", nElements);

			allData_array = printOrderParameter (dumpAtoms, dumpfile, datafile, bonds, inputVectors, currentTimestep, nElements);
			
			// for (int i = 0; i < nElements; ++i)
			// {
			// 	printf("%f %f\n", allData_array[i].distance, allData_array[i].orderParameter);
			// 	fflush (stdout);
			// 	// sleep (1);
			// }

			isTimestep = 0;
		}

		if (strstr (lineString, "ITEM: TIMESTEP"))
		{
			isTimestep = 1;
			currentDumpstep++;
			currentLine = 1;
			isNElementsSet = 0;
		}

		if (currentLine > 9 && currentLine < (9 + dumpfile.nAtoms))
		{	
			sscanf (lineString, "%d %d %f %f %f\n",
				&dumpAtoms[currentLine - 10].id,
				&dumpAtoms[currentLine - 10].atomType,
				&dumpAtoms[currentLine - 10].x,
				&dumpAtoms[currentLine - 10].y,
				&dumpAtoms[currentLine - 10].z);
		}

		// if (currentLine == (9 + dumpfile.nAtoms))
		// 	isTimestep = 1;

		currentLine++;
	}
}

int main(int argc, char const *argv[])
{
	system ("mkdir logs");

	FILE *inputDumpFile, *inputDataFile, *inputConfigFile;
	char *inputDumpFilename, *inputDataFilename, *inputConfigFilename;

	printf("%s\n", "Looking for LAMMPS trajectory file...");
	inputDumpFilename = getInputFileName ();
	printf("%s\n", "Looking for LAMMPS data file...");
	inputDataFilename = getInputFileName ();
	printf("%s\n", "Looking for input config file...");
	inputConfigFilename = getInputFileName ();

	inputDumpFile = fopen (inputDumpFilename, "r");
	inputDataFile = fopen (inputDataFilename, "r");
	inputConfigFile = fopen (inputConfigFilename, "r");

	DATA_ATOMS *atoms;
	DATA_BONDS *bonds;
	DATA_ANGLES *angles;
	DATA_DIHEDRALS *dihedrals;
	DATA_IMPROPERS *impropers;

	DATAFILE_INFO datafile;
	datafile = readData (inputDataFile, &atoms, &bonds, &angles, &dihedrals, &impropers);

	DUMPFILE_INFO dumpfile;
	dumpfile = getDumpFileInfo (inputDumpFile);

	CONFIG *inputVectors;
	inputVectors = readConfig (inputConfigFile);

	// float maxDist, hyp1, xDist = (dumpfile.xhi - dumpfile.xlo), yDist = (dumpfile.yhi - dumpfile.ylo), zDist = (dumpfile.zhi - dumpfile.zlo);
	// hyp1 = sqrt ((xDist * xDist) + (zDist * zDist));
	// maxDist = sqrt ((hyp1 * hyp1) + (yDist * yDist));

	// float binSize_dist = 1, binSize_OOP = 0.01, binSize_deg = 1;
	// int nBins_dist = (((int) maxDist) / (int) binSize_dist) + 1, nBins_OOP = (int) ((1 + 0.5) / binSize_OOP) + 1, nBins_deg = (180 / (int) binSize_deg) + 1;

	// // [degrees][distance] and [oop][distance]
	// int *distribution_degrees, *distribution_OOP, size_degrees = nBins_dist * nBins_deg, size_oop = nBins_dist * nBins_OOP;
	// distribution_OOP = (int *) calloc (size_oop, sizeof (int));
	// distribution_degrees = (int *) calloc (size_degrees, sizeof (int));

	computeOrderParameter (inputDumpFile, datafile, bonds, inputVectors);

	return 0;
}