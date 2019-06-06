#include <fstream>
#include <iomanip>
#include <random>
#include <iostream>
#include <time.h>
#include <vector>

void print_message()
{
  std::cout << "\n";
  std::cout << "#############################################################\n";
  std::cout << "\n";
  std::cout << "  Code Name: LATTICE2DPD v3\n";
  std::cout << "  Author(s): Yidong Xia (Idaho National Lab)\n";
  std::cout << "  Functions: Generate initial DPD bead locations from voxels\n";
  std::cout << "\n";
  std::cout << "#############################################################\n";
  std::cout << "\n";
}

int main(int argc, char **argv)
{
  print_message();

  clock_t timeStart = clock();

  // ========================
  // Parameter initialization
  // ========================

  // by default we use additive filling

  bool packAdd = true;

  bool outSolid = false;
  bool outFluid = false;

  unsigned int flowDirection = 0;

  // So far we consider at most two types of fluids

  unsigned int typeVoid = 0;
  unsigned int typePore = 1;
  unsigned int typeWall = 3;
  unsigned int typeWallKeep = 333;

  unsigned long ix,iy,iz;

  unsigned int cropXlo = 0, cropXhi = 100;
  unsigned int cropYlo = 0, cropYhi = 100;
  unsigned int cropZlo = 0, cropZhi = 100;

  const unsigned long nMaxTypes = 12;
  unsigned long nTypes = 0;
  unsigned long nSolidsPerLattice = 0;
  unsigned long nFluidsPerLattice = 0;
  unsigned long nMaxSolidsPerLattice = 8;
  unsigned long nMaxFluidsPerLattice = 8;

  unsigned long numSolidBeads = 0;
  unsigned long numFluidBeads = 0;

  unsigned long numSolidVox = 0;
  unsigned long numFluidVox = 0;

  unsigned long numLmpBoxBeads = 0;
  unsigned long xloLmpBox = 0;
  unsigned long xhiLmpBox = 0;
  unsigned long yloLmpBox = 0;
  unsigned long yhiLmpBox = 0;
  unsigned long zloLmpBox = 0;
  unsigned long zhiLmpBox = 0;


  unsigned long indexSolidBead = 0;
  unsigned long indexFluidBead = 0;

  unsigned long nxSolidVoxLo = 0, nxSolidVoxHi = 0;
  unsigned long nySolidVoxLo = 0, nySolidVoxHi = 0;
  unsigned long nzSolidVoxLo = 0, nzSolidVoxHi = 0;
  unsigned long nxFluidVoxLo = 0, nxFluidVoxHi = 0;
  unsigned long nyFluidVoxLo = 0, nyFluidVoxHi = 0;
  unsigned long nzFluidVoxLo = 0, nzFluidVoxHi = 0;

  unsigned long hx = 0;
  unsigned long hy = 0;
  unsigned long hz = 0;

  double xrnd, yrnd, zrnd;
  double xs, ys, zs;
  double tmp1, tmp2;
  unsigned long itmp;

  double radius = 0.5;

  // strings

  std::string skipLine;
  std::string infoFileName = "info";
  std::string ctrlFileName = "control";
  std::string inpSolidFileName = "inp_solid.dat";
  std::string inpFluidFileName = "inp_fluid.dat";
  std::string inpLAMMPSDatName = "LAMMPS_box.dat";
  std::string outSolidFileName = "out_solid.dat";
  std::string outFluidFileName = "out_fluid.dat";
  std::string outSolidDumpName = "out_solid.dump";
  std::string outFluidDumpName = "out_fluid.dump";

  // I/O files

  std::ifstream infoFile;
  std::ifstream ctrlFile;
  std::ifstream inpSolidFile;
  std::ifstream inpFluidFile;
  std::ifstream inpLAMMPSDat;
  std::ofstream outSolidFile;
  std::ofstream outFluidFile;
  std::ofstream outSolidDump;
  std::ofstream outFluidDump;

  // ==========================
  // Open and read control file
  // ==========================

  ctrlFile.open(ctrlFileName, std::ios::in);

  if (!ctrlFile.is_open())
  {
    std::cout << "\n" << "Fatal: file " << ctrlFileName << " does not exist.\n";
    std::exit(0);
  }

  // Read the flow direction indicator

  std::getline(ctrlFile,skipLine);
  ctrlFile >> flowDirection; std::getline(ctrlFile,skipLine);

  if (flowDirection > 3)
  {
    std::cout
      << "\n" << "Fatal: in file " << ctrlFileName
      << " flowDirection = " << flowDirection << " is not a valid value.\n";
    std::exit(0);
  }

  // Read the lower and higher bounds for region cropping

  std::getline(ctrlFile,skipLine);
  ctrlFile >> cropXlo >> cropXhi; std::getline(ctrlFile,skipLine);
  std::getline(ctrlFile,skipLine);
  ctrlFile >> cropYlo >> cropYhi; std::getline(ctrlFile,skipLine);
  std::getline(ctrlFile,skipLine);
  ctrlFile >> cropZlo >> cropZhi; std::getline(ctrlFile,skipLine);

  if (cropXlo < 0 || cropXhi > 100)
  {
    std::cout
      << "\n" << "Fatal: in file " << ctrlFileName
      << " cropXlo = " << cropXlo << ", "
      << " cropXhi = " << cropXhi << "\n";
    std::exit(0);
  }
  if (cropYlo < 0 || cropYhi > 100)
  {
    std::cout
      << "\n" << "Fatal: in file " << ctrlFileName
      << " cropYlo = " << cropYlo << ", "
      << " cropYhi = " << cropYhi << "\n";
    std::exit(0);
  }
  if (cropZlo < 0 || cropZhi > 100)
  {
    std::cout
      << "\n" << "Fatal: in file " << ctrlFileName
      << " cropZlo = " << cropZlo << ", "
      << " cropZhi = " << cropZhi << "\n";
    std::exit(0);
  }

  // Read the scaling coefficient in each direction

  std::getline(ctrlFile,skipLine);
  ctrlFile >> hx >> hy >> hz; std::getline(ctrlFile,skipLine);

  // Read the solid and fluid number densities

  std::getline(ctrlFile,skipLine);
  ctrlFile >> nSolidsPerLattice >> nFluidsPerLattice; std::getline(ctrlFile,skipLine);

  if(nSolidsPerLattice == 0 || nFluidsPerLattice == 0)
  {
    packAdd = false;
    std::cout << "\n" << "Packing mode: delete from full box\n\n";

    // We still need to reset them to 1 for counting pore and wall voxels
    nSolidsPerLattice = 1;
    nFluidsPerLattice = 1;
  }

  if (nSolidsPerLattice > nMaxSolidsPerLattice)
  {
    std::cout
      << "\n" << "Fatal: in file " << ctrlFileName
      << " nSolidsPerLattice = " << nSolidsPerLattice << ", "
      << " nMaxSolidsPerLattice = " << nMaxSolidsPerLattice << "\n";
    std::exit(0);
  }
  if (nFluidsPerLattice > nMaxFluidsPerLattice)
  {
    std::cout
      << "\n" << "Fatal: in file " << ctrlFileName
      << " nFluidsPerLattice = " << nFluidsPerLattice << ", "
      << " nMaxFluidsPerLattice = " << nMaxFluidsPerLattice << "\n";
    std::exit(0);
  }

  // Read the number of types of beads

  std::getline(ctrlFile,skipLine);
  ctrlFile >> nTypes; std::getline(ctrlFile,skipLine);

  if (nTypes > nMaxTypes)
  {
    std::cout
      << "\n" << "Fatal: in file " << ctrlFileName << "\n"
      << "  nTypes = " << nTypes << "\n"
      << "  is larger than\n"
      << "  nMaxTypes = " << nMaxTypes << "\n";
    std::exit(0);
  }

  // Read whether to output solid and fluid

  std::getline(ctrlFile,skipLine);
  ctrlFile >> outSolid >> outFluid; std::getline(ctrlFile,skipLine);

  ctrlFile.close();

  // =======================
  // Open and read info file
  // =======================

  // See if the file exists

  infoFile.open(infoFileName, std::ios::in);

  if (!infoFile.is_open())
  {
    std::cout
      << "\n" << "Fatal: file " << infoFileName << " does not exist.\n";
    std::exit(0);
  }

  // Get the number of solid voxels, and low & high ends of indices

  std::getline(infoFile,skipLine);
  infoFile
    >> numSolidVox;
  std::getline(infoFile,skipLine);

  std::getline(infoFile,skipLine);
  infoFile
    >> nxSolidVoxLo >> nxSolidVoxHi
    >> nySolidVoxLo >> nySolidVoxHi
    >> nzSolidVoxLo >> nzSolidVoxHi;
  std::getline(infoFile,skipLine);

  // Get the number of fluid voxels, and low & high ends of indices

  std::getline(infoFile,skipLine);
  infoFile
    >> numFluidVox;
  std::getline(infoFile,skipLine);

  std::getline(infoFile,skipLine);
  infoFile
    >> nxFluidVoxLo >> nxFluidVoxHi
    >> nyFluidVoxLo >> nyFluidVoxHi
    >> nzFluidVoxLo >> nzFluidVoxHi;
  std::getline(infoFile,skipLine);

  infoFile.close();

  // ========================================================================
  // A full 3D matrix for storing the region info is created to help identify
  // inner sides of wall voxels adjacent to pore space
  // ... void = 0 (default)
  // ... wall = 1
  // ... pore = 2
  // ========================================================================

  unsigned long nxVox = nxSolidVoxHi - nxSolidVoxLo + 1;
  unsigned long nyVox = nySolidVoxHi - nySolidVoxLo + 1;
  unsigned long nzVox = nzSolidVoxHi - nzSolidVoxLo + 1;

  unsigned long nxLattice = nxVox * hx;
  unsigned long nyLattice = nyVox * hy;
  unsigned long nzLattice = nzVox * hz;

  std::vector<std::vector<std::vector<unsigned int> > > LATTICE;
  LATTICE.resize(nxLattice+2);
  for (unsigned long i = 0; i < (nxLattice+2); i++)
  {
    LATTICE[i].resize(nyLattice+2);
    for (unsigned long j = 0; j < (nyLattice+2); j++)
      LATTICE[i][j].resize(nzLattice+2);
  }

  // ===========================
  // Display on screen and check
  // ===========================

  unsigned long xlo = nxVox * cropXlo / 100; xlo = xlo * hx;
  unsigned long xhi = nxVox * cropXhi / 100; xhi = xhi * hx;

  unsigned long ylo = nyVox * cropYlo / 100; ylo = ylo * hy;
  unsigned long yhi = nyVox * cropYhi / 100; yhi = yhi * hy;

  unsigned long zlo = nzVox * cropZlo / 100; zlo = zlo * hz;
  unsigned long zhi = nzVox * cropZhi / 100; zhi = zhi * hz;

  unsigned long numSolidLattice = numSolidVox * hx * hy * hz;
  unsigned long numFluidLattice = numFluidVox * hx * hy * hz;

  std::cout
    << std::left << std::setfill('.')
    << std::setw(40) << "numSolidVox  "       << "  " << numSolidVox  << "\n"
    << std::setw(40) << "numFluidVox  "       << "  " << numFluidVox  << "\n"
    << std::setw(40) << "nxSolidVoxLo  "      << "  " << nxSolidVoxLo << "\n"
    << std::setw(40) << "nxSolidVoxHi  "      << "  " << nxSolidVoxHi << "\n"
    << std::setw(40) << "nySolidVoxLo  "      << "  " << nySolidVoxLo << "\n"
    << std::setw(40) << "nySolidVoxHi  "      << "  " << nySolidVoxHi << "\n"
    << std::setw(40) << "nzSolidVoxLo  "      << "  " << nzSolidVoxLo << "\n"
    << std::setw(40) << "nzSolidVoxHi  "      << "  " << nzSolidVoxHi << "\n"
    << std::setw(40) << "nxFluidVoxLo  "      << "  " << nxFluidVoxLo << "\n"
    << std::setw(40) << "nxFluidVoxHi  "      << "  " << nxFluidVoxHi << "\n"
    << std::setw(40) << "nyFluidVoxLo  "      << "  " << nyFluidVoxLo << "\n"
    << std::setw(40) << "nyFluidVoxHi  "      << "  " << nyFluidVoxHi << "\n"
    << std::setw(40) << "nzFluidVoxLo  "      << "  " << nzFluidVoxLo << "\n"
    << std::setw(40) << "nzFluidVoxHi  "      << "  " << nzFluidVoxHi << "\n"
    << std::setw(40) << "nxVox  "             << "  " << nxVox        << "\n"
    << std::setw(40) << "nyVox  "             << "  " << nyVox        << "\n"
    << std::setw(40) << "nzVox  "             << "  " << nzVox        << "\n"
    << std::setw(40) << "flowDirection  "     << "  " << flowDirection<< "\n"
    << std::setw(40) << "hx  "                << "  " << hx           << "\n"
    << std::setw(40) << "hy  "                << "  " << hy           << "\n"
    << std::setw(40) << "hz  "                << "  " << hz           << "\n"
    << std::setw(40) << "xlo  "               << "  " << xlo          << "\n"
    << std::setw(40) << "xhi  "               << "  " << xhi          << "\n"
    << std::setw(40) << "ylo  "               << "  " << ylo          << "\n"
    << std::setw(40) << "yhi  "               << "  " << yhi          << "\n"
    << std::setw(40) << "zlo  "               << "  " << zlo          << "\n"
    << std::setw(40) << "zhi  "               << "  " << zhi          << "\n"
    << std::setw(40) << "nTypes  "            << "  " << nTypes       << "\n"
    << std::setw(40) << "nSolidsPerLattice  " << "  " << nSolidsPerLattice<< "\n"
    << std::setw(40) << "nFluidsPerLattice  " << "  " << nFluidsPerLattice<< "\n"
    << std::setw(40) << "outSolid  "          << "  " << (outSolid ? "Yes" : "No") << "\n"
    << std::setw(40) << "outFluid  "          << "  " << (outFluid ? "Yes" : "No") << "\n"
    << "\n";

  // =====================================
  // Open and read input solid voxel files
  // =====================================

  inpSolidFile.open(inpSolidFileName, std::ios::in);

  if (!inpSolidFile.is_open())
  {
    std::cout
      << "\n" << "Fatal: file " << inpSolidFileName << " does not exist.\n";
    std::exit(0);
  }
  else
    std::cout << "\n" << "Reading file " << inpSolidFileName << " ... ";

  for (unsigned long i = 0; i < numSolidVox; i++)
  {
    inpSolidFile >> ix >> iy >> iz;
    std::getline(inpSolidFile,skipLine);

    bool isOutOfRange = false;

    if (ix < nxSolidVoxLo || ix > nxSolidVoxHi) isOutOfRange = true;
    if (iy < nySolidVoxLo || iy > nySolidVoxHi) isOutOfRange = true;
    if (iz < nzSolidVoxLo || iz > nzSolidVoxHi) isOutOfRange = true;
    if (isOutOfRange)
    {
      std::cout << "\n\n"
        << "Error: in file " << inpSolidFileName << ", Line " << (i+1) << ", "
        << "(ix,iy,iz) = (" << ix << "," << iy << "," << iz << ") out of range!\n";
      std::exit(0);
    }

    // Remember each index starts from ONE
    unsigned long I = ix - nxSolidVoxLo + 1;
    unsigned long J = iy - nySolidVoxLo + 1;
    unsigned long K = iz - nzSolidVoxLo + 1;

    unsigned long I0 = (I-1) * hx + 1;
    unsigned long I1 =     I * hx;
    unsigned long J0 = (J-1) * hy + 1;
    unsigned long J1 =     J * hy;
    unsigned long K0 = (K-1) * hz + 1;
    unsigned long K1 =     K * hz;

    for (unsigned long ii = I0; ii <= I1; ii++)
      for (unsigned long jj = J0; jj <= J1; jj++)
        for (unsigned long kk = K0; kk <= K1; kk++)
          LATTICE[ii][jj][kk] = typeWall;
  }

  inpSolidFile.close();
  std::cout << "complete\n ";

  // =====================================
  // Open and read input fluid voxel files
  // =====================================

  inpFluidFile.open(inpFluidFileName, std::ios::in);

  if (!inpFluidFile.is_open())
  {
    std::cout
      << "\n" << "Fatal: file " << inpFluidFileName << " does not exist.\n";
    std::exit(0);
  }
  else
    std::cout << "\n" << "Reading file " << inpFluidFileName << " ... ";

  for (unsigned long i = 0; i < numFluidVox; i++)
  {
    inpFluidFile >> ix >> iy >> iz;
    std::getline(inpFluidFile,skipLine);

    bool isOutOfRange = false;

    if (ix < nxFluidVoxLo || ix > nxFluidVoxHi) isOutOfRange = true;
    if (iy < nyFluidVoxLo || iy > nyFluidVoxHi) isOutOfRange = true;
    if (iz < nzFluidVoxLo || iz > nzFluidVoxHi) isOutOfRange = true;
    if (isOutOfRange)
    {
      std::cout << "\n\n"
        << "Error: in file " << inpFluidFileName << ", Line " << (i+1) << ", "
        << "(ix,iy,iz) = (" << ix << "," << iy << "," << iz << ") out of range!\n";
      std::exit(0);
    }

    // Remember each index starts from ONE
    unsigned long I = ix - nxSolidVoxLo + 1;
    unsigned long J = iy - nySolidVoxLo + 1;
    unsigned long K = iz - nzSolidVoxLo + 1;

    unsigned long I0 = (I-1) * hx + 1;
    unsigned long I1 =     I * hx;
    unsigned long J0 = (J-1) * hy + 1;
    unsigned long J1 =     J * hy;
    unsigned long K0 = (K-1) * hz + 1;
    unsigned long K1 =     K * hz;

    for (unsigned long ii = I0; ii <= I1; ii++)
      for (unsigned long jj = J0; jj <= J1; jj++)
        for (unsigned long kk = K0; kk <= K1; kk++)
          LATTICE[ii][jj][kk] = typePore;
  }

  inpFluidFile.close();
  std::cout << "complete\n";

  // =========================================================================
  // Find the total numbers of solid and fluid particles in the cropped region
  // =========================================================================

  for (unsigned long i = (xlo+1); i <= xhi; i++)
    for (unsigned long j = (ylo+1); j <= yhi; j++)
      for (unsigned long k = (zlo+1); k <= zhi; k++)
      {
        if ((i == (xlo+1)) || (i == xhi) || (j == (ylo+1)) || (j == yhi) || (k == (zlo+1)) || (k == zhi))
        {
          if (LATTICE[i][j][k] == typeWall)
          {
            LATTICE[i][j][k] = typeWallKeep;
            numSolidBeads += nSolidsPerLattice;
          }
          else if (LATTICE[i][j][k] == typePore)
          {
            LATTICE[i][j][k] = typeWall;
            numSolidBeads += nSolidsPerLattice;
          }
        }
        else
        {
          if (LATTICE[i][j][k] == typeWall)
          {
            bool fillSolids = false;

            for (unsigned long ii = i-1; ii <= (i+1); ii++)
            {
              for (unsigned long jj = j-1; jj <= (j+1); jj++)
              {
                for (unsigned long kk = k-1; kk <= (k+1); kk++)
                {
                  if(LATTICE[ii][jj][kk] == typePore)
                  {
                    fillSolids = true;
                    numSolidBeads += nSolidsPerLattice;
                    break;
                  }
                }
                if (fillSolids) break;
              }
              if (fillSolids) break;
            }

            if (!fillSolids)
              LATTICE[i][j][k] = typeVoid;
          }
          else if (LATTICE[i][j][k] == typePore)
            numFluidBeads += nFluidsPerLattice;
        }
      }

  if (flowDirection == 0)
  {
  }
  else if (flowDirection == 1)
  {
    for (unsigned long j = (ylo+1); j <= yhi; j++)
      for (unsigned long k = (zlo+1); k <= zhi; k++)
      {
        unsigned long i;

        i = xlo+1;

        if ( (LATTICE[i  ][j  ][k  ] == typeWall) &&
             (LATTICE[i+1][j-1][k-1] != typeVoid) &&
             (LATTICE[i+1][j-1][k  ] != typeVoid) &&
             (LATTICE[i+1][j-1][k+1] != typeVoid) &&
             (LATTICE[i+1][j  ][k-1] != typeVoid) &&
             (LATTICE[i+1][j  ][k  ] != typeVoid) &&
             (LATTICE[i+1][j  ][k+1] != typeVoid) &&
             (LATTICE[i+1][j+1][k-1] != typeVoid) &&
             (LATTICE[i+1][j+1][k  ] != typeVoid) &&
             (LATTICE[i+1][j+1][k+1] != typeVoid) )
        {
          if ( (j > (ylo+1)) && (j < yhi) && (k > (zlo+1)) && (k < zhi) )
          {
            LATTICE[i][j][k] = typePore;
            numFluidBeads += nFluidsPerLattice;
            numSolidBeads -= nSolidsPerLattice;
          }
        }
        else if (LATTICE[i][j][k] == typeWallKeep)
        {
          LATTICE[i][j][k] = typeWall;
        }
        else if (LATTICE[i][j][k] == typeVoid)
        {
          LATTICE[i][j][k] = typeWall;
          numSolidBeads += nSolidsPerLattice;
        }

        i = xhi;

        if ( (LATTICE[i  ][j  ][k  ] == typeWall) &&
             (LATTICE[i-1][j-1][k-1] != typeVoid) &&
             (LATTICE[i-1][j-1][k  ] != typeVoid) &&
             (LATTICE[i-1][j-1][k+1] != typeVoid) &&
             (LATTICE[i-1][j  ][k-1] != typeVoid) &&
             (LATTICE[i-1][j  ][k  ] != typeVoid) &&
             (LATTICE[i-1][j  ][k+1] != typeVoid) &&
             (LATTICE[i-1][j+1][k-1] != typeVoid) &&
             (LATTICE[i-1][j+1][k  ] != typeVoid) &&
             (LATTICE[i-1][j+1][k+1] != typeVoid) )
        {
          if ( (j > (ylo+1)) && (j < yhi) && (k > (zlo+1)) && (k < zhi) )
          {
            LATTICE[i][j][k] = typePore;
            numFluidBeads += nFluidsPerLattice;
            numSolidBeads -= nSolidsPerLattice;
          }
        }
        else if (LATTICE[i][j][k] == typeWallKeep)
        {
          LATTICE[i][j][k] = typeWall;
        }
        else if (LATTICE[i][j][k] == typeVoid)
        {
          LATTICE[i][j][k] = typeWall;
          numSolidBeads += nSolidsPerLattice;
        }
      }
  }
  else if (flowDirection == 2)
  {
    for (unsigned long i = (xlo+1); i <= xhi; i++)
      for (unsigned long k = (zlo+1); k <= zhi; k++)
      {
        unsigned long j;

        j = ylo+1;

        if ( (LATTICE[i  ][j  ][k  ] == typeWall) &&
             (LATTICE[i-1][j+1][k-1] != typeVoid) &&
             (LATTICE[i-1][j+1][k  ] != typeVoid) &&
             (LATTICE[i-1][j+1][k+1] != typeVoid) &&
             (LATTICE[i  ][j+1][k-1] != typeVoid) &&
             (LATTICE[i  ][j+1][k  ] != typeVoid) &&
             (LATTICE[i  ][j+1][k+1] != typeVoid) &&
             (LATTICE[i+1][j+1][k-1] != typeVoid) &&
             (LATTICE[i+1][j+1][k  ] != typeVoid) &&
             (LATTICE[i+1][j+1][k+1] != typeVoid) )
        {
          if ( (i > (xlo+1)) && (i < xhi) && (k > (zlo+1)) && (k < zhi) )
          {
            LATTICE[i][j][k] = typePore;
            numFluidBeads += nFluidsPerLattice;
            numSolidBeads -= nSolidsPerLattice;
          }
        }
        else if (LATTICE[i][j][k] == typeWallKeep)
        {
          LATTICE[i][j][k] = typeWall;
        }
        else if (LATTICE[i][j][k] == typeVoid)
        {
          LATTICE[i][j][k] = typeWall;
          numSolidBeads += nSolidsPerLattice;
        }

        j = yhi;

        if ( (LATTICE[i  ][j  ][k  ] == typeWall) &&
             (LATTICE[i-1][j-1][k-1] != typeVoid) &&
             (LATTICE[i-1][j-1][k  ] != typeVoid) &&
             (LATTICE[i-1][j-1][k+1] != typeVoid) &&
             (LATTICE[i  ][j-1][k-1] != typeVoid) &&
             (LATTICE[i  ][j-1][k  ] != typeVoid) &&
             (LATTICE[i  ][j-1][k+1] != typeVoid) &&
             (LATTICE[i+1][j-1][k-1] != typeVoid) &&
             (LATTICE[i+1][j-1][k  ] != typeVoid) &&
             (LATTICE[i+1][j-1][k+1] != typeVoid) )
        {
          if ( (i > (xlo+1)) && (i < xhi) && (k > (zlo+1)) && (k < zhi) )
          {
            LATTICE[i][j][k] = typePore;
            numFluidBeads += nFluidsPerLattice;
            numSolidBeads -= nSolidsPerLattice;
          }
        }
        else if (LATTICE[i][j][k] == typeWallKeep)
        {
          LATTICE[i][j][k] = typeWall;
        }
        else if (LATTICE[i][j][k] == typeVoid)
        {
          LATTICE[i][j][k] = typeWall;
          numSolidBeads += nSolidsPerLattice;
        }
      }
  }
  else if (flowDirection == 3)
  {
    for (unsigned long i = (xlo+1); i <= xhi; i++)
      for (unsigned long j = (ylo+1); j <= yhi; j++)
      {
        unsigned long k;

        k = zlo+1;

        if ( (LATTICE[i  ][j  ][k  ] == typeWall) &&
             (LATTICE[i-1][j-1][k+1] != typeVoid) &&
             (LATTICE[i-1][j  ][k+1] != typeVoid) &&
             (LATTICE[i-1][j+1][k+1] != typeVoid) &&
             (LATTICE[i  ][j-1][k+1] != typeVoid) &&
             (LATTICE[i  ][j  ][k+1] != typeVoid) &&
             (LATTICE[i  ][j+1][k+1] != typeVoid) &&
             (LATTICE[i+1][j-1][k+1] != typeVoid) &&
             (LATTICE[i+1][j  ][k+1] != typeVoid) &&
             (LATTICE[i+1][j+1][k+1] != typeVoid) )
        {
          if ( (i > (xlo+1)) && (i < xhi) && (j > (ylo+1)) && (j < yhi) )
          {
            LATTICE[i][j][k] = typePore;
            numFluidBeads += nFluidsPerLattice;
            numSolidBeads -= nSolidsPerLattice;
          }
        }
        else if (LATTICE[i][j][k] == typeWallKeep)
        {
          LATTICE[i][j][k] = typeWall;
        }
        else if (LATTICE[i][j][k] == typeVoid)
        {
          LATTICE[i][j][k] = typeWall;
          numSolidBeads += nSolidsPerLattice;
        }

        k = zhi;

        if ( (LATTICE[i  ][j  ][k  ] == typeWall) &&
             (LATTICE[i-1][j-1][k-1] != typeVoid) &&
             (LATTICE[i-1][j  ][k-1] != typeVoid) &&
             (LATTICE[i-1][j+1][k-1] != typeVoid) &&
             (LATTICE[i  ][j-1][k-1] != typeVoid) &&
             (LATTICE[i  ][j  ][k-1] != typeVoid) &&
             (LATTICE[i  ][j+1][k-1] != typeVoid) &&
             (LATTICE[i+1][j-1][k-1] != typeVoid) &&
             (LATTICE[i+1][j  ][k-1] != typeVoid) &&
             (LATTICE[i+1][j+1][k-1] != typeVoid) )
        {
          if ( (i > (xlo+1)) && (i < xhi) && (j > (ylo+1)) && (j < yhi) )
          {
            LATTICE[i][j][k] = typePore;
            numFluidBeads += nFluidsPerLattice;
            numSolidBeads -= nSolidsPerLattice;
          }
        }
        else if (LATTICE[i][j][k] == typeWallKeep)
        {
          LATTICE[i][j][k] = typeWall;
        }
        else if (LATTICE[i][j][k] == typeVoid)
        {
          LATTICE[i][j][k] = typeWall;
          numSolidBeads += nSolidsPerLattice;
        }
      }
  }

  numSolidVox = numSolidLattice / hx / hy / hz;
  numFluidVox = numFluidLattice / hx / hy / hz;

  // ========================================================================
  // Get the number of solid and fluid beads in the case of packing by delete
  // ========================================================================

  std::vector<std::vector<double> > BEADS_XYZ;
  std::vector<unsigned long> BEADS_TYPE;
  std::vector<bool> BEADS_KEEP;

  if (!packAdd)
  {
    inpLAMMPSDat.open(inpLAMMPSDatName, std::ios::in);

    if (!inpLAMMPSDat.is_open())
    {
      std::cout
        << "\n" << "Fatal: file " << inpLAMMPSDatName << " does not exist.\n";
      std::exit(0);
    }
    else
      std::cout << "\n" << "Reading file " << inpLAMMPSDatName << "\n";

    std::getline(inpLAMMPSDat,skipLine);
    std::getline(inpLAMMPSDat,skipLine);

    inpLAMMPSDat >> numLmpBoxBeads;
    std::getline(inpLAMMPSDat,skipLine);

    std::getline(inpLAMMPSDat,skipLine);

    inpLAMMPSDat >> tmp1 >> tmp2;
    std::getline(inpLAMMPSDat,skipLine);
    tmp1 = tmp1 + 0.5 - (tmp1<0); xloLmpBox = (int)tmp1;
    tmp2 = tmp2 + 0.5 - (tmp2<0); xhiLmpBox = (int)tmp2;

    inpLAMMPSDat >> tmp1 >> tmp2;
    std::getline(inpLAMMPSDat,skipLine);
    tmp1 = tmp1 + 0.5 - (tmp1<0); yloLmpBox = (int)tmp1;
    tmp2 = tmp2 + 0.5 - (tmp2<0); yhiLmpBox = (int)tmp2;

    inpLAMMPSDat >> tmp1 >> tmp2;
    std::getline(inpLAMMPSDat,skipLine);
    tmp1 = tmp1 + 0.5 - (tmp1<0); zloLmpBox = (int)tmp1;
    tmp2 = tmp2 + 0.5 - (tmp2<0); zhiLmpBox = (int)tmp2;

    std::cout
      << std::left << std::setfill('.')
      << std::setw(40) << "numLmpBoxBeads  "       << "  " << numLmpBoxBeads  << "\n"
      << std::setw(40) << "xloLmpBox  "            << "  " << xloLmpBox       << "\n"
      << std::setw(40) << "xhiLmpBox  "            << "  " << xhiLmpBox       << "\n"
      << std::setw(40) << "yloLmpBox  "            << "  " << yloLmpBox       << "\n"
      << std::setw(40) << "yhiLmpBox  "            << "  " << yhiLmpBox       << "\n"
      << std::setw(40) << "zloLmpBox  "            << "  " << zloLmpBox       << "\n"
      << std::setw(40) << "zhiLmpBox  "            << "  " << zhiLmpBox       << "\n";

    if (xhiLmpBox != (xhi-xlo))
    {
      std::cout << "\n\n" << "Error: xhiLmpBox is not equal to (xhi - xlo) ... Please check!\n";
      std::exit(0);
    }
    if (yhiLmpBox != (yhi-ylo))
    {
      std::cout << "\n\n" << "Error: yhiLmpBox is not equal to (yhi - ylo) ... Please check!\n";
      std::exit(0);
    }
    if (zhiLmpBox != (zhi-zlo))
    {
      std::cout << "\n\n" << "Error: zhiLmpBox is not equal to (zhi - zlo) ... Please check!\n";
      std::exit(0);
    }

    for (unsigned int i = 0; i < 11; i++)
      std::getline(inpLAMMPSDat,skipLine);

    BEADS_KEEP.resize(numLmpBoxBeads, false);
    BEADS_TYPE.resize(numLmpBoxBeads);
    BEADS_XYZ.resize(numLmpBoxBeads);

    unsigned long j = 0;
    unsigned long jmax = 0;

    for (unsigned long i = 0; i < numLmpBoxBeads; i++)
    {
      inpLAMMPSDat >> j;
      if (j > jmax) jmax = j;

      BEADS_XYZ[j-1].resize(3);

      inpLAMMPSDat >> BEADS_TYPE[j-1] >> BEADS_XYZ[j-1][0] >> BEADS_XYZ[j-1][1] >> BEADS_XYZ[j-1][2];
      std::getline(inpLAMMPSDat,skipLine);

      if (((i+1)%1000) == 0)
        std::cout << "\r " << std::fixed << std::setw(5) << std::setprecision(3) << (100.0*(double)(i+1)/(double)(numLmpBoxBeads)) << " % complete" << std::flush;
    }
    if (jmax != numLmpBoxBeads)
    {
      std::cout << "\n" << "Fatal: max particle index in input LAMMPS data is " << jmax << " ... It is NOT equal to numLmpBoxBeads.\n";
      std::exit(0);
    }

    inpLAMMPSDat.close();

    // Reset these numbers

    numSolidBeads = 0;
    numFluidBeads = 0;

    // Count the number of fluid and solid beads to be output

    std::cout << "\n\n" << "Counting the total number of solid and fluid beads for output:\n";

    unsigned long int_xp;
    unsigned long int_yp;
    unsigned long int_zp;

    for (unsigned long l = 0; l < numLmpBoxBeads; l++)
    {
      int_xp = (unsigned long)(BEADS_XYZ[l][0]);
      int_yp = (unsigned long)(BEADS_XYZ[l][1]);
      int_zp = (unsigned long)(BEADS_XYZ[l][2]);

      unsigned long i = xlo + 1 + int_xp;
      unsigned long j = ylo + 1 + int_yp;
      unsigned long k = zlo + 1 + int_zp;

      if (i > xhi)
      {
        std::cout << "\n\n" << "Error: lattice index i is greater than xhi ... Please check!\n";
        std::cout << "i   = " << i   << "\n";
        std::cout << "xhi = " << xhi << "\n";
        std::exit(0);
      }
      if (j > yhi)
      {
        std::cout << "\n\n" << "Error: lattice index j is greater than yhi ... Please check!\n";
        std::cout << "j   = " << j   << "\n";
        std::cout << "yhi = " << yhi << "\n";
        std::exit(0);
      }
      if (k > zhi)
      {
        std::cout << "\n\n" << "Error: lattice index k is greater than zhi ... Please check!\n";
        std::cout << "k   = " << k   << "\n";
        std::cout << "zhi = " << zhi << "\n";
        std::exit(0);
      }

      if (LATTICE[i][j][k] != typeVoid)
      {
        BEADS_KEEP[l] = true;

        if ((LATTICE[i][j][k] == typeWall) || (LATTICE[i][j][k] == typeWallKeep))
        {
          numSolidBeads ++;
          BEADS_TYPE[l] = typeWall;
        }
        else if ((LATTICE[i][j][k] == typePore))
        {
          numFluidBeads ++;
          BEADS_TYPE[l] = typePore;
        }
        else
        {
          std::cout << "\n\n" << "Error: exceptional lattice type ... Please check!\n";
          std::cout << "i = " << i << "\n";
          std::cout << "j = " << j << "\n";
          std::cout << "k = " << k << "\n";
          std::cout << "lattice type = " << LATTICE[i][j][k] << "\n";
          std::exit(0);
        }
      }
      if (((l+1)%1000) == 0)
        std::cout << "\r " << std::fixed << std::setw(5) << std::setprecision(3) << (100.0*(double)(l+1)/(double)(numLmpBoxBeads)) << " % complete" << std::flush;
    }
  }

  // ===========================
  // Output the LAMMPS data file
  // ===========================

  std::cout
    << "\n"
    << "\n" << "Writing output files: "
    << "\n" << "  " << outSolidFileName
    << "\n" << "  " << outSolidDumpName
    << "\n" << "  " << outFluidFileName
    << "\n" << "  " << outFluidDumpName
    << "\n";

  outSolidFile.open(outSolidFileName, std::ios::out);
  outSolidFile.precision(16);
  outSolidFile << std::scientific;

  outSolidFile
    << "LAMMPS data file for read_data\n"
    << "\n"
    << numSolidBeads << " atoms\n"
    << nTypes << " atom types\n"
    << "\n"
    << 0 << " " << (xhi-xlo) << " xlo xhi\n"
    << 0 << " " << (yhi-ylo) << " ylo yhi\n"
    << 0 << " " << (zhi-zlo) << " zlo zhi\n"
    << "\n"
    << "Atoms\n"
    << "\n";

  outFluidFile.open(outFluidFileName, std::ios::out);
  outFluidFile.precision(16);
  outFluidFile << std::scientific;

  outFluidFile
    << "LAMMPS data file for read_data\n"
    << "\n"
    << numFluidBeads << " atoms\n"
    << nTypes << " atom types\n"
    << "\n"
    << 0 << " " << (xhi-xlo) << " xlo xhi\n"
    << 0 << " " << (yhi-ylo) << " ylo yhi\n"
    << 0 << " " << (zhi-zlo) << " zlo zhi\n"
    << "\n"
    << "Atoms\n"
    << "\n";

  // LAMMPS dump file format

  outSolidDump.open(outSolidDumpName, std::ios::out);
  outSolidDump.precision(16);
  outSolidDump << std::scientific;

  outSolidDump
    << "ITEM: TIMESTEP\n"
    << "0\n"
    << "ITEM: NUMBER OF ATOMS\n"
    << numSolidBeads << "\n"
    << "ITEM: BOX BOUNDS ff ff ff\n"
    << 0 << " " << (xhi-xlo) << "\n"
    << 0 << " " << (yhi-ylo) << "\n"
    << 0 << " " << (zhi-zlo) << "\n"
    << "ITEM: ATOMS id type x y z radius\n";

  outFluidDump.open(outFluidDumpName, std::ios::out);
  outFluidDump.precision(16);
  outFluidDump << std::scientific;

  outFluidDump
    << "ITEM: TIMESTEP\n"
    << "0\n"
    << "ITEM: NUMBER OF ATOMS\n"
    << numFluidBeads << "\n"
    << "ITEM: BOX BOUNDS ff ff ff\n"
    << 0 << " " << (xhi-xlo) << "\n"
    << 0 << " " << (yhi-ylo) << "\n"
    << 0 << " " << (zhi-zlo) << "\n"
    << "ITEM: ATOMS id type x y z radius\n";

  // ==============================================================================
  // Choose the way to do particle packing and write particle coordinates into file
  // ==============================================================================

  if (packAdd)
  {
    // Open the DPD thermodynamically-equilibrated unit cube distribution

    std::string iNameSolidsRhoUnitDist = "LAMMPS_rho" + std::to_string(nSolidsPerLattice) + "_unit_dist.dump";
    std::ifstream iFileSolidsRhoUnitDist;
    iFileSolidsRhoUnitDist.open(iNameSolidsRhoUnitDist, std::ios::in);
    if (nSolidsPerLattice > 1)
    {
      if (!iFileSolidsRhoUnitDist.is_open())
      {
        std::cout << "\n" << "Fatal: file " << iNameSolidsRhoUnitDist << " does not exist.\n";
        std::exit(0);
      }
      else
        std::cout << "Open file: " << iNameSolidsRhoUnitDist << std::endl;
    }

    std::string iNameFluidsRhoUnitDist = "LAMMPS_rho" + std::to_string(nFluidsPerLattice) + "_unit_dist.dump";
    std::ifstream iFileFluidsRhoUnitDist;
    iFileFluidsRhoUnitDist.open(iNameFluidsRhoUnitDist, std::ios::in);
    if (nFluidsPerLattice > 1)
    {
      if (!iFileFluidsRhoUnitDist.is_open())
      {
        std::cout << "\n" << "Fatal: file " << iNameFluidsRhoUnitDist << " does not exist.\n";
        std::exit(0);
      }
      else
        std::cout << "Open file: " << iNameFluidsRhoUnitDist << std::endl;
    }

    for (unsigned long i = (xlo+1); i <= xhi; i++)
      for (unsigned long j = (ylo+1); j <= yhi; j++)
        for (unsigned long k = (zlo+1); k <= zhi; k++)
        {
          if ((LATTICE[i][j][k] == typeWall) || (LATTICE[i][j][k] == typeWallKeep))
          {
            if (nSolidsPerLattice == 1)
            {
              xrnd = 0.5 + double(i-1-xlo);
              yrnd = 0.5 + double(j-1-ylo);
              zrnd = 0.5 + double(k-1-zlo);

              indexSolidBead ++;

              outSolidFile << indexSolidBead << " " << typeWall << " " << xrnd << " " << yrnd << " " << zrnd << "\n";
              outSolidDump << indexSolidBead << " " << typeWall << " " << xrnd << " " << yrnd << " " << zrnd << " " << radius << "\n";
            }
            else
            {
              // At the end of file, close file and open again
              if (iFileSolidsRhoUnitDist.eof())
              {
                iFileSolidsRhoUnitDist.close();
                iFileSolidsRhoUnitDist.open(iNameSolidsRhoUnitDist, std::ios::in);
              }

              // Skip 9 lines
              std::getline(iFileSolidsRhoUnitDist,skipLine);
              std::getline(iFileSolidsRhoUnitDist,skipLine);
              std::getline(iFileSolidsRhoUnitDist,skipLine);
              std::getline(iFileSolidsRhoUnitDist,skipLine);
              std::getline(iFileSolidsRhoUnitDist,skipLine);
              std::getline(iFileSolidsRhoUnitDist,skipLine);
              std::getline(iFileSolidsRhoUnitDist,skipLine);
              std::getline(iFileSolidsRhoUnitDist,skipLine);
              std::getline(iFileSolidsRhoUnitDist,skipLine);

              for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
              {
                iFileSolidsRhoUnitDist >> itmp >> itmp >> xs >> ys >> zs;
                std::getline(iFileSolidsRhoUnitDist,skipLine);

                // take care of periodic boundary condition
                if (xs < 0.001) xs = 0.999 - (0.001 - xs);
                if (xs > 0.999) xs = 0.001 + (xs - 0.999);
                if (ys < 0.001) ys = 0.999 - (0.001 - ys);
                if (ys > 0.999) ys = 0.001 + (ys - 0.999);
                if (zs < 0.001) zs = 0.999 - (0.001 - zs);
                if (zs > 0.999) zs = 0.001 + (zs - 0.999);

                xrnd = xs + double(i-1-xlo);
                yrnd = ys + double(j-1-ylo);
                zrnd = zs + double(k-1-zlo);

                indexSolidBead ++;

                outSolidFile << indexSolidBead << " " << typeWall << " " << xrnd << " " << yrnd << " " << zrnd << "\n";
                outSolidDump << indexSolidBead << " " << typeWall << " " << xrnd << " " << yrnd << " " << zrnd << " " << radius << "\n";
              }
            }
          }
          else if ((LATTICE[i][j][k] == typePore) && outFluid)
          {
            if (nFluidsPerLattice == 1)
            {
              xrnd = 0.5 + double(i-1-xlo);
              yrnd = 0.5 + double(j-1-ylo);
              zrnd = 0.5 + double(k-1-zlo);

              indexFluidBead ++;

              outFluidFile << indexFluidBead << " " << typePore << " " << xrnd << " " << yrnd << " " << zrnd << "\n";
              outFluidDump << indexFluidBead << " " << typePore << " " << xrnd << " " << yrnd << " " << zrnd << " " << radius << "\n";
            }
            else
            {
              // At the end of file, close file and open again
              if (iFileFluidsRhoUnitDist.eof())
              {
                iFileFluidsRhoUnitDist.close();
                iFileFluidsRhoUnitDist.open(iNameFluidsRhoUnitDist, std::ios::in);
              }

              // Skip 9 lines
              std::getline(iFileFluidsRhoUnitDist,skipLine);
              std::getline(iFileFluidsRhoUnitDist,skipLine);
              std::getline(iFileFluidsRhoUnitDist,skipLine);
              std::getline(iFileFluidsRhoUnitDist,skipLine);
              std::getline(iFileFluidsRhoUnitDist,skipLine);
              std::getline(iFileFluidsRhoUnitDist,skipLine);
              std::getline(iFileFluidsRhoUnitDist,skipLine);
              std::getline(iFileFluidsRhoUnitDist,skipLine);
              std::getline(iFileFluidsRhoUnitDist,skipLine);

              for (unsigned int irnd = 0; irnd < nFluidsPerLattice; irnd++)
              {
                iFileFluidsRhoUnitDist >> itmp >> itmp >> xs >> ys >> zs;
                std::getline(iFileFluidsRhoUnitDist,skipLine);

                // take care of periodic boundary condition
                if (xs < 0.001) xs = 0.999 - (0.001 - xs);
                if (xs > 0.999) xs = 0.001 + (xs - 0.999);
                if (ys < 0.001) ys = 0.999 - (0.001 - ys);
                if (ys > 0.999) ys = 0.001 + (ys - 0.999);
                if (zs < 0.001) zs = 0.999 - (0.001 - zs);
                if (zs > 0.999) zs = 0.001 + (zs - 0.999);

                xrnd = xs + double(i-1-xlo);
                yrnd = ys + double(j-1-ylo);
                zrnd = zs + double(k-1-zlo);

                indexFluidBead ++;

                outFluidFile << indexFluidBead << " " << typePore << " " << xrnd << " " << yrnd << " " << zrnd << "\n";
                outFluidDump << indexFluidBead << " " << typePore << " " << xrnd << " " << yrnd << " " << zrnd << " " << radius << "\n";
              }
            }
          }
        }

    if (iFileSolidsRhoUnitDist.is_open())
      iFileSolidsRhoUnitDist.close();

    if (iFileFluidsRhoUnitDist.is_open())
      iFileFluidsRhoUnitDist.close();
  }
  else
  {
    for (unsigned long l = 0; l < numLmpBoxBeads; l++)
    {
      if (BEADS_KEEP[l])
      {
        if((BEADS_TYPE[l] == typeWall || BEADS_TYPE[l] == typeWallKeep) && outSolid)
        {
          indexSolidBead ++;

          outSolidFile << indexSolidBead << " " << typeWall << " " << BEADS_XYZ[l][0] << " " << BEADS_XYZ[l][1] << " " << BEADS_XYZ[l][2] << "\n";
          outSolidDump << indexSolidBead << " " << typeWall << " " << BEADS_XYZ[l][0] << " " << BEADS_XYZ[l][1] << " " << BEADS_XYZ[l][2] << " " << radius << "\n";
        }
        else if (BEADS_TYPE[l] == typePore && outFluid)
        {
          indexFluidBead ++;

          outFluidFile << indexFluidBead << " " << typePore << " " << BEADS_XYZ[l][0] << " " << BEADS_XYZ[l][1] << " " << BEADS_XYZ[l][2] << "\n";
          outFluidDump << indexFluidBead << " " << typePore << " " << BEADS_XYZ[l][0] << " " << BEADS_XYZ[l][1] << " " << BEADS_XYZ[l][2] << " " << radius << "\n";
        }
      }
      if (((l+1)%1000) == 0)
        std::cout << "\r " << std::fixed << std::setw(5) << std::setprecision(3) << (100.0*(double)(l+1)/(double)(numLmpBoxBeads)) << " % complete" << std::flush;
    }
  }

  // ======================
  // Close active I/O files
  // ======================

  outSolidFile.close();
  outFluidFile.close();
  outSolidDump.close();
  outFluidDump.close();

  // =======================
  // Final data sanity check
  // =======================

  std::cout << "\n\n" << "In the cropped region:\n\n";
  std::cout
    << std::left << std::setfill('.')
    << std::setw(40) << "xlo  "            << "  " << "0"      << "\n"
    << std::setw(40) << "xhi  "            << "  " << xhi-xlo  << "\n"
    << std::setw(40) << "ylo  "            << "  " << "0"      << "\n"
    << std::setw(40) << "yhi  "            << "  " << yhi-ylo  << "\n"
    << std::setw(40) << "zlo  "            << "  " << "0"      << "\n"
    << std::setw(40) << "zhi  "            << "  " << zhi-zlo  << "\n"
    << std::setw(40) << "numSolidBeads  "  << "  " << numSolidBeads  << "\n"
    << std::setw(40) << "indexSolidBead  " << "  " << indexSolidBead << "\n"
    << std::setw(40) << "numFluidBeads  "  << "  " << numFluidBeads  << "\n"
    << std::setw(40) << "indexFluidBead  " << "  " << indexFluidBead << "\n";

  // ====================
  // Report the wall time
  // ====================

  double timeElapsed = (double)(clock() - timeStart) / CLOCKS_PER_SEC;
  std::cout << "\n" << "Time elapsed in seconds: " << timeElapsed << "\n";
}
