/*
	Copyright Markus Kühbach, 2015-2017

	TopologyTracer is an MPI-parallel datamining tool for the conducting of spatiotemporal 
	analyses of microstructural element dynamics. Its purpose is the quantification of correlations
	--- spatial and temporal --- between individual microstructural elements and their 
	higher-order neighboring elements. Specifically, its current functionalities allow to 
	study the volume evolution of individual grains over time and set their growth history into 
	relation to the evolution of the higher-order neighbors. The tool is unique insofar as it 
	allows processing these individual surveys in a parallelized manner. Thus, enabling the 
	post-processing of so far intractable large datasets.

	The source code was developed by Markus Kühbach during his PhD time with Luis A. Barrales-Mora 
	and Günter Gottstein at the Institute of Physical Metallurgy and Metal Physics with RWTH Aachen University. 
	Being now with the Max-Planck-Institut fur Eisenforschung GmbH in Dusseldorf, I maintain the code, 
	though at disregular intervals. Nonetheless, feel free to utilize the tool, do not hesitate contacting 
	me for sharing thoughts, suggesting improvements, or reporting your experiences.
	markus.kuehbach at rwth-aachen.de and m.kuehbach at mpie.de


	The authors gratefully acknowledge the financial support from the Deutsche Forschungsgemeinschaft
	(DFG) within the Reinhart Koselleck-Project (GO 335/44-1) and computing time grants kindly provided
	by RWTH Aachen University and the FZ Jülich within the scope of the JARAHPC project JARA0076.


	This file is part of TopologyTracer.

	TopologyTracer is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	TopologyTracer is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with SCORE.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "TopologyTracer2D3D_Settings.h"

//initialization of the setting variables
unsigned int Settings::SimID = 0;
double Settings::LocalDatabaseMaximumSize = 0.0;
unsigned int Settings::SnapshotFirst = 0;
unsigned int Settings::SnapshotOffset = 0;
unsigned int Settings::SnapshotLast = 0;
unsigned int Settings::Dimensionality = 2;
double Settings::MPIReadBlockLength = MPIIO_DEFAULT_BLOCKLENGTH;
double Settings::HAGBMobility = 0.0;
double Settings::HAGBEnergy = 0.0;
double Settings::DislocEnPerM = 0.0;
double Settings::PhysicalDomainSize = 0.0;
unsigned int Settings::LookupMaxGrainIDRange = EMPIRICAL_LOOKUP_MAXIDRANGE;
unsigned int Settings::LargestGrainID = EMPIRICAL_LARGEST_GRAINID;
unsigned int Settings::MemRegionsX = 1;
unsigned int Settings::MemRegionsY = 1;
unsigned int Settings::MemRegionsZ = 1;
unsigned int Settings::MaxNumberOfKShells = DEFAULT_KSHELL_MAX;
double Settings::HAGBDetectionThreshold = DEFAULT_HAGBTHRESHOLD;
double Settings::DisoriAngleBinMin = DEFAULT_MINIMUM_DISORI_FCC;
double Settings::DisoriAngleBinMax = DEFAULT_MAXIMUM_DISORI_FCC;
double Settings::DisoriAngleBinWidth = DEFAULT_DISORI_BINWIDTH;
double Settings::SeeBinMin = DEFAULT_MINIMUM_SEE;
double Settings::SeeBinMax = DEFAULT_MAXIMUM_SEE;
double Settings::SeeBinWidth = DEFAULT_SEE_BINWIDTH;
double Settings::GSDBinMin = DEFAULT_MINIMUM_GSD;
double Settings::GSDBinMax = DEFAULT_MAXIMUM_GSD;
double Settings::GSDBinWidth = DEFAULT_GSD_BINWIDTH;
bool Settings::DeveloperMode = false;
bool Settings::ProbeWorkPartMode = false;
bool Settings::ProbeBoundaryContact = false;
double Settings::LongRangeRadiusMin = 0.0;
double Settings::LongRangeRadiusIncr = 0.0;
double Settings::LongRangeRadiusMax = 0.0;
double Settings::InitialDomainEdgeLength = 0.0;

bool Settings::AnalyzeTrajectoriesForward = false;
bool Settings::AnalyzeTrajectoriesBackward = false;
bool Settings::AnalyzeSizeGainVsMatrixBackward = false;
bool Settings::AnalyzeMaxSizeGainForward = false;
bool Settings::AnalyzeMeanDrvForceSEEForward = false;
bool Settings::AnalyzeKNN = false;
bool Settings::AnalyzeGrainSizeQuantiles = false;
bool Settings::AnalyzeMODF = false;
bool Settings::AnalyzeSEE = false;
bool Settings::AnalyzeGSD = false;
bool Settings::AnalyzeDrvForceSEE = false;
bool Settings::AnalyzeClassicalNucModels = false;
bool Settings::AnalyzeApproxRXFraction = false;
bool Settings::AnalyzeAbnormalGrains = false;
bool Settings::AnalyzeTopologyDifferenceForward = false;

unsigned int Settings::AnalyzeTrajectoriesForwardMode = E_FORWARD_ALL;
bool Settings::ComputeCurvatureTJP = false;
bool Settings::TranslateBinary2GNU = false;

E_ANALYSIS_MODE Settings::AnalysisMode = E_TRACKING_PARALLEL;
E_COARSENING_MODEL Settings::SimModel = E_LEVELSET;

E_CAPTRACKING_MODE Settings::CapillaryTrackingMode = E_CAP_AVERAGE;
E_CAPTRACKING_WHICH Settings::CapillaryTrackingTargets = E_CAP_ALL;


string Settings::TargetGIDFromFilename;
string Settings::KNNDataFromFilename;
string Settings::GNUDataFromFilename;
string Settings::UDSDataFromFilename;
string Settings::GrainBndContactFilename;

unsigned int Settings::get_GSDNQuantiles(void)
{
	return STATISTICS_VOL_NQUANTILES;
}

unsigned int Settings::get_SEEBinCount(void)
{
	double binend = ((Settings::SeeBinMax - Settings::SeeBinMin) / Settings::SeeBinWidth);
	unsigned int nbins = 1 + ((unsigned int) std::ceil(binend)) + 1; //+1 for all smaller than Settings::SeeBinMin and +1 for all larger than Settings::SeeBinMax
	return nbins;
}


unsigned int Settings::get_MODFBinCount(void)
{
	double binend = (Settings::DisoriAngleBinMax - Settings::DisoriAngleBinMin) / Settings::DisoriAngleBinWidth;
	unsigned int nbins = 1 + ((unsigned int) std::ceil(binend)) + 1; //+1 for all smaller smallest, +1 for all larger range end

	return nbins;
}


unsigned int Settings::get_GSDBinCount(void)
{
	double binend = (Settings::GSDBinMax - Settings::GSDBinMin) / Settings::GSDBinWidth;
	unsigned int nbins = 1 + ((unsigned int) std::ceil(binend)) + 1; //+1 for all smaller smallest, +1 for all larger range end

	return nbins;
}


double Settings::get_SEEBinEnd(void) {
	double binend = (Settings::SeeBinMax - Settings::SeeBinMin) / Settings::SeeBinWidth;
	return binend;
}


double Settings::get_MODFBinEnd(void) {
	double binend = (Settings::DisoriAngleBinMax - Settings::DisoriAngleBinMin) / Settings::DisoriAngleBinWidth;
	return binend;
}


double Settings::get_GSDBinEnd(void) {
	double binend = (Settings::GSDBinMax - Settings::GSDBinMin) / Settings::GSDBinWidth;
	return binend;
}


//definition of the functions
void Settings::readXML(string filename) {

	//find the desired .xml file
	if ( 0 == filename.compare("") )
		filename = string("TopologyTracer2D3D_Parameter.xml");
	ifstream file( filename );
	if (file.fail()) {
		throw runtime_error(string("Unable to locate file ") + filename);
	}
	stringstream contents;
	contents << file.rdbuf();
	string xmlDocument(contents.str());

	xml_document<> tree;
	tree.parse<0>(&xmlDocument[0]);

	xml_node<>* rootNode = tree.first_node();
	if (0 != strcmp(rootNode->name(), "Parameters")) {
		throw runtime_error("undefined parameters file!");
	}

	//Read all the required parameters
	/*if (0 != rootNode->first_node("SimulationID")) {
		SimID = std::stoul(
				rootNode->first_node("SimulationID")->value());
	}*/

	/*if (0 != rootNode->first_node("SimulationModel")) {
		SimModel = (E_COARSENING_MODEL) std::stoi(
				rootNode->first_node("SimulationModel")->value());

		if (SimModel >= E_LEVELSET)
			SimModel = E_LEVELSET;
	}*/
	unsigned int mode = E_UNKNOWN_MODE;
	if (0 != rootNode->first_node("AnalysisMode")) {
		mode = std::stoul( rootNode->first_node("AnalysisMode")->value());
	}

	if ( mode == E_TRACKING_PARALLEL ) 
		AnalysisMode = E_TRACKING_PARALLEL;
	else if ( mode == E_TRACKING_SEQUENTIAL )
		AnalysisMode = E_TRACKING_SEQUENTIAL;
	else if ( mode == E_KNNCORRELATION ) 
		AnalysisMode = E_KNNCORRELATION;
	else if ( mode == E_NUCSITES_LONGRANGE2D )
		AnalysisMode = E_NUCSITES_LONGRANGE2D;
	else if ( mode == E_NUCSITES_LONGRANGE3D )
		AnalysisMode = E_NUCSITES_LONGRANGE3D;
	else if ( mode == E_APPROX_CURVATURE_CONTOUR )
		AnalysisMode = E_APPROX_CURVATURE_CONTOUR;
	else if ( mode == E_CAPILLARY_ACTIVITY )
		AnalysisMode = E_CAPILLARY_ACTIVITY;
	else 
		AnalysisMode = E_UNKNOWN_MODE;

	/*switch ( mode ) //##MK::somewhat buggy
	{
		case E_TRACKING_PARALLEL :
			AnalysisMode = E_TRACKING_PARALLEL; break;
		case E_TRACKING_SEQUENTIAL :
			AnalysisMode = E_TRACKING_SEQUENTIAL; break;
		case E_KNNCORRELATION :
			AnalysisMode = E_KNNCORRELATION; break;
		case E_NUCSITES_LONGRANGE2D :
			AnalysisMode = E_NUCSITES_LONGRANGE2D; break;
		case E_NUCSITES_LONGRANGE3D :
			AnalysisMode = E_NUCSITES_LONGRANGE3D; break;
		case E_APPROX_CURVATURE_CONTOUR :
			AnalysisMode = E_APPROX_CURVATURE_CONTOUR; break;
		case E_CAPILLARY_ACTIVITY :
			AnalysisMode = E_CAPILLARY_ACTIVITY; break;
		default : 
			AnalysisMode = E_UNKNOWN_MODE; break;
	}*/

	if (0 != rootNode->first_node("TargetGrainIDs")) {
		TargetGIDFromFilename = rootNode->first_node("TargetGrainIDs")->value();
	}
	if (0 != rootNode->first_node("KNNFile")) {
		KNNDataFromFilename = rootNode->first_node("KNNFile")->value();
	}
	if (0 != rootNode->first_node("GNUFile")) {
		GNUDataFromFilename = rootNode->first_node("GNUFile")->value();
	}
	if (0 != rootNode->first_node("UDSFile")) {
		UDSDataFromFilename = rootNode->first_node("UDSFile")->value();
	}
	if (0 != rootNode->first_node("GrainBoundaryContact")) {
		GrainBndContactFilename = rootNode->first_node("GrainBoundaryContact")->value();
	}
	
	SimModel = E_LEVELSET;
	if (0 != rootNode->first_node("DeveloperMode")) {
		DeveloperMode = std::stoul(
				rootNode->first_node("DeveloperMode")->value());
	}	
	if (0 != rootNode->first_node("OnlyProbeTheWorkPartitioning")) {
		ProbeWorkPartMode = std::stoul(
				rootNode->first_node("OnlyProbeTheWorkPartitioning")->value());
	}
	if (0 != rootNode->first_node("ProbeBoundaryContact")) {
		ProbeBoundaryContact = std::stoul(
			rootNode->first_node("ProbeBoundaryContact")->value());
	}
	if (0 != rootNode->first_node("LocalDatabaseMaximumSize")) {
		LocalDatabaseMaximumSize = std::stod(
				rootNode->first_node("LocalDatabaseMaximumSize")->value());
	}

	if (0 != rootNode->first_node("SnapshotFirst")) {
		SnapshotFirst = std::stoul(
				rootNode->first_node("SnapshotFirst")->value());
	}
	if (0 != rootNode->first_node("SnapshotOffset")) {
		SnapshotOffset = std::stoul(
				rootNode->first_node("SnapshotOffset")->value());
	}
	if (0 != rootNode->first_node("SnapshotLast")) {
		SnapshotLast = std::stoul(
				rootNode->first_node("SnapshotLast")->value());
	}
	if (0 != rootNode->first_node("Dimensionality")) {
		Dimensionality = std::stoul(
				rootNode->first_node("Dimensionality")->value());
	}
	if (0 != rootNode->first_node("MPIReadBlockLength")) {
		MPIReadBlockLength = std::stod(
				rootNode->first_node("MPIReadBlockLength")->value());
	}
	if (0 != rootNode->first_node("HAGBMobility")) {
		HAGBMobility = std::stod(
			rootNode->first_node("HAGBMobility")->value());
	}
	if (0 != rootNode->first_node("HAGBEnergy")) {
		HAGBEnergy = std::stod(
				rootNode->first_node("HAGBEnergy")->value());
	}
	if (0 != rootNode->first_node("DislocEnPerM")) {
		DislocEnPerM = std::stod(
				rootNode->first_node("DislocEnPerM")->value());
	}
	if (0 != rootNode->first_node("PhysDomainEdgeLength")) {
		PhysicalDomainSize = std::stod(
				rootNode->first_node("PhysDomainEdgeLength")->value());
	}

	if (0 != rootNode->first_node("MaxIDRange")) {
		LookupMaxGrainIDRange = std::stoul(
				rootNode->first_node("MaxIDRange")->value());
	}
   	if (0 != rootNode->first_node("LargestGrainID")) {
		LargestGrainID = std::stoul(
				rootNode->first_node("LargestGrainID")->value());
	}     
	if (0 != rootNode->first_node("MemRegionsX")) {
		MemRegionsX = std::stoul(
				rootNode->first_node("MemRegionsX")->value());
	}	
	if (0 != rootNode->first_node("MemRegionsY")) {
		MemRegionsY = std::stoul(
				rootNode->first_node("MemRegionsY")->value());
	}	
	if (0 != rootNode->first_node("MemRegionsZ")) {
		MemRegionsZ = std::stoul(
				rootNode->first_node("MemRegionsZ")->value());
	}	
	if (0 != rootNode->first_node("MaximumNumberOfKShells")) {
		MaxNumberOfKShells = std::stoul(
				rootNode->first_node("MaximumNumberOfKShells")->value());
	}
	if (0 != rootNode->first_node("HAGBDetectionThreshold")) {
		HAGBDetectionThreshold = std::stod(
				rootNode->first_node("HAGBDetectionThreshold")->value());
	}

	if (0 != rootNode->first_node("DisoriAngleBinningMin")) {
		DisoriAngleBinMin = std::stod(
				rootNode->first_node("DisoriAngleBinningMin")->value());
	}
	if (0 != rootNode->first_node("DisoriAngleBinningMax")) {
		DisoriAngleBinMax = std::stod(
				rootNode->first_node("DisoriAngleBinningMax")->value());
	}
	if (0 != rootNode->first_node("DisoriAngleBinningWidth")) {
		DisoriAngleBinWidth = std::stod(
				rootNode->first_node("DisoriAngleBinningWidth")->value());
	}
	if (0 != rootNode->first_node("StoredElasticEnergyBinningMin")) {
		SeeBinMin = std::stod( rootNode->first_node("StoredElasticEnergyBinningMin")->value());
	}
	if (0 != rootNode->first_node("StoredElasticEnergyBinningMax")) {
		SeeBinMax = std::stod( rootNode->first_node("StoredElasticEnergyBinningMax")->value());
	}
	if (0 != rootNode->first_node("StoredElasticEnergyBinningWidth")) {
		SeeBinWidth = std::stod( rootNode->first_node("StoredElasticEnergyBinningWidth")->value());
	}
	if (0 != rootNode->first_node("GSDBinningMin")) {
		GSDBinMin = std::stod( rootNode->first_node("GSDBinningMin")->value());
	}
	if (0 != rootNode->first_node("GSDBinningMax")) {
		GSDBinMax = std::stod( rootNode->first_node("GSDBinningMax")->value());
	}
	if (0 != rootNode->first_node("GSDBinningWidth")) {
		GSDBinWidth = std::stod( rootNode->first_node("GSDBinningWidth")->value());
	}

	if (0 != rootNode->first_node("PersistenceRadiusMin")) {
		LongRangeRadiusMin = std::stod(
				rootNode->first_node("PersistenceRadiusMin")->value());
	}
	if (0 != rootNode->first_node("PersistenceRadiusIncr")) {
		LongRangeRadiusIncr = std::stod(
				rootNode->first_node("PersistenceRadiusIncr")->value());
	}
	if (0 != rootNode->first_node("PersistenceRadiusMax")) {
		LongRangeRadiusMax = std::stod(
				rootNode->first_node("PersistenceRadiusMax")->value());
	}
	if (0 != rootNode->first_node("InitialDomainEdgeLength")) {
		InitialDomainEdgeLength = std::stod(
			rootNode->first_node("InitialDomainEdgeLength")->value());
	}
	

	//which analyses to perform
	if (0 != rootNode->first_node("AnalyzeTrajectoriesForward")) {
		AnalyzeTrajectoriesForward = std::stoul(
				rootNode->first_node("AnalyzeTrajectoriesForward")->value());
	}
	if (0 != rootNode->first_node("AnalyzeTrajectoriesBackward")) {
		AnalyzeTrajectoriesBackward = std::stoul(
				rootNode->first_node("AnalyzeTrajectoriesBackward")->value());
	}
	if (0 != rootNode->first_node("AnalyzeSizeGainVsMatrixBackward")) {
		AnalyzeSizeGainVsMatrixBackward = std::stoul(
				rootNode->first_node("AnalyzeSizeGainVsMatrixBackward")->value());
	}
	if (0 != rootNode->first_node("AnalyzeMaxSizeGainForward")) {
		AnalyzeMaxSizeGainForward = std::stoul(
			rootNode->first_node("AnalyzeMaxSizeGainForward")->value());
	}
	if (0 != rootNode->first_node("AnalyzeMeanDrivingForceSEEForward")) {
		AnalyzeMeanDrvForceSEEForward = std::stoul(
			rootNode->first_node("AnalyzeMeanDrivingForceSEEForward")->value());
	}
	if (0 != rootNode->first_node("AnalyzeKNN")) {
		AnalyzeKNN = std::stoul(
				rootNode->first_node("AnalyzeKNN")->value());
	}
	if (0 != rootNode->first_node("AnalyzeGrainSizeQuantiles")) {
		AnalyzeGrainSizeQuantiles = std::stoul(
				rootNode->first_node("AnalyzeGrainSizeQuantiles")->value());
	}
	if (0 != rootNode->first_node("AnalyzeMODF")) {
		AnalyzeMODF = std::stoul(
			rootNode->first_node("AnalyzeMODF")->value());
	}
	if(0 != rootNode->first_node("AnalyzeSEE")) {
		AnalyzeSEE = std::stoul(
			rootNode->first_node("AnalyzeSEE")->value());
	}
	if(0 != rootNode->first_node("AnalyzeGSD")) {
		AnalyzeGSD = std::stoul(
			rootNode->first_node("AnalyzeGSD")->value());
	}
	if (0 != rootNode->first_node("AnalyzeDrivingForceSEE")) {
		AnalyzeDrvForceSEE = std::stoul(
			rootNode->first_node("AnalyzeDrivingForceSEE")->value());
	}
	if (0 != rootNode->first_node("AnalyzeClassicalNucModels")) {
		AnalyzeClassicalNucModels = std::stoul(
				rootNode->first_node("AnalyzeClassicalNucModels")->value());
	}
	if (0 != rootNode->first_node("AnalyzeApproxRXFraction")) {
		AnalyzeApproxRXFraction = std::stoul(
				rootNode->first_node("AnalyzeApproxRXFraction")->value());
	}
	if (0 != rootNode->first_node("AnalyzeAbnormalGrains")) {
		AnalyzeAbnormalGrains = std::stoul(
				rootNode->first_node("AnalyzeAbnormalGrains")->value());
	}
	if (0 != rootNode->first_node("AnalyzeTopologyDifferenceForward")) {
		AnalyzeTopologyDifferenceForward = std::stoul(
				rootNode->first_node("AnalyzeTopologyDifferenceForward")->value());
	}
	if (0 != rootNode->first_node("AnalyzeTrajectoriesForwardMode")) {
		AnalyzeTrajectoriesForwardMode = std::stoul(
				rootNode->first_node("AnalyzeTrajectoriesForwardMode")->value());
	}

	if (0 != rootNode->first_node("ComputeCurvatureAlsoAtTJP")) {
		ComputeCurvatureTJP = std::stoul(
				rootNode->first_node("ComputeCurvatureAlsoAtTJP")->value());
	}
	if (0 != rootNode->first_node("TranslateBinary2GNU")) {
		TranslateBinary2GNU = std::stoul(
				rootNode->first_node("TranslateBinary2GNU")->value());
	}

	mode = E_CAP_ALL;
	if (0 != rootNode->first_node("CapillaryActivityTargets")) {
		mode = std::stoul( rootNode->first_node("CapillaryActivityTargets")->value() );
	}
	if ( mode == E_CAP_SEL ) 
		CapillaryTrackingTargets = E_CAP_SEL;
	else
		CapillaryTrackingTargets = E_CAP_ALL;

	mode = E_CAP_AVERAGE;
	if (0 != rootNode->first_node("CapillaryActivityMode")) {
		mode = std::stoul( rootNode->first_node("CapillaryActivityMode")->value() );
	}
	if ( mode == E_CAP_HISTORY )
		CapillaryTrackingMode = E_CAP_HISTORY;
	else 
		CapillaryTrackingMode = E_CAP_AVERAGE;


	if (0 != rootNode->first_node("TranslateBinary2GNU")) {
		TranslateBinary2GNU = std::stoul(
			rootNode->first_node("TranslateBinary2GNU")->value());
	}
	
}

bool Settings::checkUserInput( void )
{
	/*checks the user input for validity and good sense*/
	//cout << setprecision(18) << Settings::LocalDatabaseMaximumSize << "\t\t" << Settings::SnapshotFirst << ";" <<  Settings::SnapshotOffset << ";" << Settings::SnapshotLast << endl;

	if ( Settings::LocalDatabaseMaximumSize <= 0.0 || Settings::LocalDatabaseMaximumSize > LOCALDBMAXIMUMSIZE ) return false;
	if ( Settings::SnapshotFirst > Settings::SnapshotLast || Settings::SnapshotOffset < 1 ) return false;
	//if ( Settings::SnapshotFirst % Settings::SnapshotOffset != 0 || Settings::SnapshotLast % Settings::SnapshotOffset != 0  ) return false;
	if ( Settings::Dimensionality < 2 || Settings::Dimensionality > 3 ) return false;
	if ( Settings::MPIReadBlockLength <= 0.0 ) return false;
	//if ( Settings::HAGBMobility <= EPSILON ) return false;
	if ( Settings::HAGBEnergy <= EPSILON ) return false;
	if ( Settings::DislocEnPerM <= EPSILON ) return false;
	//if ( Settings::PhysicalDomainSize <= EPSILON ) return false;
	if ( Settings::LargestGrainID >= EMPIRICAL_LARGEST_GRAINID ) return false;
	if ( Settings::LookupMaxGrainIDRange < EMPIRICAL_LOOKUP_MINIDRANGE ) return false;
	if ( Settings::MaxNumberOfKShells > DEFAULT_KSHELL_MAX ) return false;
	if ( Settings::HAGBDetectionThreshold < DEFAULT_HAGBTHRESHOLD_MIN || Settings:: HAGBDetectionThreshold > DEFAULT_HAGBTHRESHOLD_MAX )
		HAGBDetectionThreshold = DEFAULT_HAGBTHRESHOLD;

	if ( Settings::DisoriAngleBinMin <= 0.0 || Settings::DisoriAngleBinMin > DEFAULT_MAXIMUM_DISORI_FCC )
		Settings::DisoriAngleBinMin = DEFAULT_MINIMUM_DISORI_FCC;
	if ( Settings::DisoriAngleBinMax <= 0.0 || Settings::DisoriAngleBinMax > DEFAULT_MAXIMUM_DISORI_FCC )
		Settings::DisoriAngleBinMax = DEFAULT_MAXIMUM_DISORI_FCC;
	if ( Settings::DisoriAngleBinWidth <= 0.0 || Settings::DisoriAngleBinWidth > 2.0 ) 
		DisoriAngleBinWidth = DEFAULT_DISORI_BINWIDTH;
	if (get_MODFBinCount() < 2)
		return false;

	if ( Settings::SeeBinMin <= 0.0 || Settings::SeeBinMin > DEFAULT_MINIMUM_SEE ) 
		Settings::SeeBinMin = DEFAULT_MINIMUM_SEE;
	if ( Settings::SeeBinMax <= 0.0 || Settings::SeeBinMax > DEFAULT_MAXIMUM_SEE ) 
		Settings::SeeBinMax = DEFAULT_MAXIMUM_SEE;
	if ( Settings::SeeBinWidth <= 0.0 || Settings::SeeBinWidth < Settings::SeeBinMin || Settings::SeeBinWidth > Settings::SeeBinMax ) 
		Settings::SeeBinWidth = DEFAULT_SEE_BINWIDTH;
	if ( get_SEEBinCount() < 2)
		return false;

	if ( Settings::GSDBinMin <= 0.0 || Settings::GSDBinMin > DEFAULT_MINIMUM_GSD ) 
		Settings::GSDBinMin = DEFAULT_MINIMUM_GSD;
	if ( Settings::GSDBinMax <= 0.0 || Settings::GSDBinMax > DEFAULT_MAXIMUM_GSD ) 
		Settings::GSDBinMax = DEFAULT_MAXIMUM_GSD;
	if ( Settings::GSDBinWidth <= 0.0 || Settings::GSDBinWidth < Settings::GSDBinMin || Settings::GSDBinWidth > Settings::GSDBinMax ) 
		Settings::GSDBinWidth = DEFAULT_GSD_BINWIDTH;
	if ( get_GSDBinCount() < 2)
		return false;


	if ( Settings::AnalysisMode == E_NUCSITES_LONGRANGE2D || Settings::AnalysisMode == E_NUCSITES_LONGRANGE3D ) {
		if ( Settings::LongRangeRadiusMin <= 0.0 || Settings::LongRangeRadiusMin >= 1.0 )	return false;
		if ( Settings::LongRangeRadiusIncr <= 0.0 || Settings::LongRangeRadiusIncr >= 1.0 )	return false;
		if ( Settings::LongRangeRadiusMax <= 0.0 || Settings::LongRangeRadiusMax >= 1.0 )	return false;
		if ( Settings::LongRangeRadiusMin > Settings::LongRangeRadiusMax ) return false;
		if ( Settings::InitialDomainEdgeLength < EPSILON ) return false;
	}

	//check for invalid analysis modes
	if ( Settings::AnalyzeTrajectoriesForward == true ) {
		if ( Settings::AnalyzeTrajectoriesForwardMode != E_FORWARD_ALL && Settings::AnalyzeTrajectoriesForwardMode != E_FORWARD_SELECTED ) return false;
	}

	//deg2rad
	HAGBDetectionThreshold = DEG2RAD(HAGBDetectionThreshold);
	DisoriAngleBinMin = DEG2RAD(DisoriAngleBinMin);
	DisoriAngleBinMax = DEG2RAD(DisoriAngleBinMax);
	DisoriAngleBinWidth = DEG2RAD(DisoriAngleBinWidth);

	return true;
}



