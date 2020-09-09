#include "GeaneEigenStorageUtils.hh"

namespace gm2strawtracker {

  void GeaneEigenStorageUtils::PrepareEigenForStorage(std::unique_ptr<DecayVertexArtRecordCollection>& recordsPtr) {

    for (auto & record : *recordsPtr) {
      
      //Eigen::MatrixXd objects
      for (int i(0); i<record.covMatrix.size(); i++) record.covMatrixData.push_back(record.covMatrix(i));
      record.covMatrix.resize(0,0);

      for (auto step : record.steps) {
	for (int i(0); i<step.jacobian.size(); i++) step.jacobianData.push_back(step.jacobian(i));
	step.jacobian.resize(0,0);
      }

    }

  }

  //-----------------------------------------------------------------------------------------------------------
  // Prepare a collection of GEANE objects with data members converted to eigen data members
  //-----------------------------------------------------------------------------------------------------------
  void GeaneEigenStorageUtils::PrepareEigenForStorage(std::unique_ptr<TrackArtRecordCollection>& recordsPtr)
  {

    // Loop over all records in collection
    for(auto & record : *recordsPtr){

      // Eigen::MatrixXd objects
      for(int i = 0; i < record.geaneHits.covarianceTotalInverse.size(); i++) record.geaneHits.covarianceTotalInverseData.push_back(record.geaneHits.covarianceTotalInverse(i));
      record.geaneHits.covarianceTotalInverse.resize(0,0);

      for(auto & matrix : record.geaneHits.geaneErrorMatrices){
        std::vector<double> data;
        for(int i = 0; i < matrix.size(); i++) data.push_back(matrix(i));
        record.geaneHits.geaneErrorMatricesData.push_back(data);
      }
      record.geaneHits.geaneErrorMatrices.clear();

      for(auto & matrix : record.geaneHits.extendedTransportMatrixBegToEnd){
        std::vector<double> data;
        for(int i = 0; i < matrix.size(); i++) data.push_back(matrix(i));
        record.geaneHits.extendedTransportMatrixBegToEndData.push_back(data);
      }
      record.geaneHits.extendedTransportMatrixBegToEnd.clear();

    }

    return;
  }

 
  //-----------------------------------------------------------------------------------------------------------
  // Prepare a collection of GEANE objects with data members converted to eigen data members
  //-----------------------------------------------------------------------------------------------------------
  void GeaneEigenStorageUtils::PrepareEigenForStorage(std::unique_ptr<TrackDetailArtRecordCollection>& recordsPtr)
  {

    // Loop over all records in collection
    for(auto & record : *recordsPtr){

      // Eigen::MatrixXd objects
      for(int i = 0; i < record.geaneHits.covarianceTotalInverse.size(); i++) record.geaneHits.covarianceTotalInverseData.push_back(record.geaneHits.covarianceTotalInverse(i));
      record.geaneHits.covarianceTotalInverse.resize(0,0);

      for(int i = 0; i < record.geaneHits.extendedCombinedTransportMatricesTranspose.size(); i++) record.geaneHits.extendedCombinedTransportMatricesTransposeData.push_back(record.geaneHits.extendedCombinedTransportMatricesTranspose(i));
      record.geaneHits.extendedCombinedTransportMatricesTranspose.resize(0,0);

      for(int i = 0; i < record.geaneHits.extendedReducedMatrix.size(); i++) record.geaneHits.extendedReducedMatrixData.push_back(record.geaneHits.extendedReducedMatrix(i));
      record.geaneHits.extendedReducedMatrix.resize(0,0);

      for(int i = 0; i < record.geaneHits.extendedReducedMatrixInverse.size(); i++) record.geaneHits.extendedReducedMatrixInverseData.push_back(record.geaneHits.extendedReducedMatrixInverse(i));
      record.geaneHits.extendedReducedMatrixInverse.resize(0,0);

      for(int i = 0; i < record.geaneHits.extendedModifiedReducedMatrix.size(); i++) record.geaneHits.extendedModifiedReducedMatrixData.push_back(record.geaneHits.extendedModifiedReducedMatrix(i));
      record.geaneHits.extendedModifiedReducedMatrix.resize(0,0);

      // std::vector<Eigen::MatrixXd>
      for(auto & matrix : record.geaneHits.geaneTransportMatrices){
        std::vector<double> data;
        for(int i = 0; i < matrix.size(); i++) data.push_back(matrix(i));
        record.geaneHits.geaneTransportMatricesData.push_back(data);
      }
      record.geaneHits.geaneTransportMatrices.clear();

      for(auto & matrix : record.geaneHits.geaneErrorMatrices){
        std::vector<double> data;
        for(int i = 0; i < matrix.size(); i++) data.push_back(matrix(i));
        record.geaneHits.geaneErrorMatricesData.push_back(data);
      }
      record.geaneHits.geaneErrorMatrices.clear();

      for(auto & matrix : record.geaneHits.extendedTransportMatrixBegToEnd){
        std::vector<double> data;
        for(int i = 0; i < matrix.size(); i++) data.push_back(matrix(i));
        record.geaneHits.extendedTransportMatrixBegToEndData.push_back(data);
      }
      record.geaneHits.extendedTransportMatrixBegToEnd.clear();

      // std::vector<Eigen::VectorXd>
      for(auto & vector : record.geaneHits.paramPredictedInUVEigen){
        std::vector<double> data;
        for(int i = 0; i < vector.size(); i++) data.push_back(vector(i));
        record.geaneHits.paramPredictedInUVEigenData.push_back(data);
      }
      record.geaneHits.paramPredictedInUVEigen.clear();
    }

    return;
  }

 
  //-----------------------------------------------------------------------------------------------------------
  // Update a collection of GEANE objects with data members converted to eigen data members
  //-----------------------------------------------------------------------------------------------------------
  TrackArtRecordCollection GeaneEigenStorageUtils::ReadEigenFromStorage(const TrackArtRecordCollection& records)
  {
    // Declare collection that we're going to return (can't update in analyser as requires const argument)
    TrackArtRecordCollection readCollection;

    // Loop over all records in collection
    for(auto & record : records){
      TrackArtRecord readRecord = record;
      TrackArtRecord updatedReadRecord = ReadEigenFromStorage(readRecord);          
      readCollection.push_back(updatedReadRecord);
    }

    // Return the updated collection
    return readCollection;
  }
 
  TrackDetailArtRecordCollection GeaneEigenStorageUtils::ReadEigenFromStorage(const TrackDetailArtRecordCollection& records)
  {
    // Declare collection that we're going to return (can't update in analyser as requires const argument)
    TrackDetailArtRecordCollection readCollection;

    // Loop over all records in collection
    for(auto & record : records){
      TrackDetailArtRecord readRecord = record;
      TrackDetailArtRecord updatedReadRecord = ReadEigenFromStorage(readRecord);          
      readCollection.push_back(updatedReadRecord);
    }

    // Return the updated collection
    return readCollection;
  }


  //-----------------------------------------------------------------------------------------------------------
  // Update a GEANE object with data members converted to eigen data members
  //-----------------------------------------------------------------------------------------------------------
  TrackDetailArtRecord GeaneEigenStorageUtils::ReadEigenFromStorage(const TrackDetailArtRecord& record)
  {
    // copy
    auto readRecord = record;

    // Eigen::MatrixXd
    int covarianceTotalInverseRows = sqrt(readRecord.geaneHits.covarianceTotalInverseData.size());
    int covarianceTotalInverseCols = sqrt(readRecord.geaneHits.covarianceTotalInverseData.size());
    readRecord.geaneHits.covarianceTotalInverse = Eigen::MatrixXd(covarianceTotalInverseRows,covarianceTotalInverseCols);
    for(unsigned int i = 0; i < readRecord.geaneHits.covarianceTotalInverseData.size(); i++){
       readRecord.geaneHits.covarianceTotalInverse(i) = readRecord.geaneHits.covarianceTotalInverseData.at(i);
    }
    readRecord.geaneHits.covarianceTotalInverseData.clear();

    // this object is size Nx5, so the matrix needs to be filled differently
    int extendedCombinedTransportMatricesTransposeRows = 5;
    int extendedCombinedTransportMatricesTransposeCols = readRecord.trackNumPlanesHit;
    readRecord.geaneHits.extendedCombinedTransportMatricesTranspose = Eigen::MatrixXd(extendedCombinedTransportMatricesTransposeRows,extendedCombinedTransportMatricesTransposeCols);
    for(unsigned int i = 0; i < readRecord.geaneHits.extendedCombinedTransportMatricesTransposeData.size(); i++){
      readRecord.geaneHits.extendedCombinedTransportMatricesTranspose(i) = readRecord.geaneHits.extendedCombinedTransportMatricesTransposeData.at(i);
    }
    readRecord.geaneHits.extendedCombinedTransportMatricesTransposeData.clear();

    int extendedReducedMatrixRows = sqrt(readRecord.geaneHits.extendedReducedMatrixData.size());
    int extendedReducedMatrixCols = sqrt(readRecord.geaneHits.extendedReducedMatrixData.size());
    readRecord.geaneHits.extendedReducedMatrix = Eigen::MatrixXd(extendedReducedMatrixRows,extendedReducedMatrixCols);
    for(unsigned int i = 0; i < readRecord.geaneHits.extendedReducedMatrixData.size(); i++){
      readRecord.geaneHits.extendedReducedMatrix(i) = readRecord.geaneHits.extendedReducedMatrixData.at(i);
    }
    readRecord.geaneHits.extendedReducedMatrixData.clear();

    int extendedReducedMatrixInverseRows = sqrt(readRecord.geaneHits.extendedReducedMatrixInverseData.size());
    int extendedReducedMatrixInverseCols = sqrt(readRecord.geaneHits.extendedReducedMatrixInverseData.size());
    readRecord.geaneHits.extendedReducedMatrixInverse = Eigen::MatrixXd(extendedReducedMatrixInverseRows,extendedReducedMatrixInverseCols);
    for(unsigned int i = 0; i < readRecord.geaneHits.extendedReducedMatrixInverseData.size(); i++){
      readRecord.geaneHits.extendedReducedMatrixInverse(i) = readRecord.geaneHits.extendedReducedMatrixInverseData.at(i);
    }
    readRecord.geaneHits.extendedReducedMatrixInverseData.clear();

    int extendedModifiedReducedMatrixRows = sqrt(readRecord.geaneHits.extendedModifiedReducedMatrixData.size());
    int extendedModifiedReducedMatrixCols = sqrt(readRecord.geaneHits.extendedModifiedReducedMatrixData.size());
    readRecord.geaneHits.extendedModifiedReducedMatrix = Eigen::MatrixXd(extendedModifiedReducedMatrixRows,extendedModifiedReducedMatrixCols);
    for(unsigned int i = 0; i < readRecord.geaneHits.extendedModifiedReducedMatrixData.size(); i++){
      readRecord.geaneHits.extendedModifiedReducedMatrix(i) = readRecord.geaneHits.extendedModifiedReducedMatrixData.at(i);
    }            
    readRecord.geaneHits.extendedModifiedReducedMatrixData.clear();

    // vector<Eigen::MatrixXd>
    for(auto & matrixData : record.geaneHits.geaneTransportMatricesData){
      int matrixRows = sqrt(matrixData.size());
      int matrixCols = sqrt(matrixData.size());
      Eigen::MatrixXd	matrix(matrixRows,matrixCols);
      for(unsigned int i = 0; i < matrixData.size(); i++){
       	matrix(i) = matrixData.at(i);
      }
      readRecord.geaneHits.geaneTransportMatrices.push_back(matrix);
    }
    readRecord.geaneHits.geaneTransportMatricesData.clear();

    for(auto & matrixData : record.geaneHits.geaneErrorMatricesData){
      int matrixRows = sqrt(matrixData.size());
      int matrixCols = sqrt(matrixData.size());
      Eigen::MatrixXd matrix(matrixRows,matrixCols);
      for(unsigned int i = 0; i < matrixData.size(); i++){
        matrix(i) = matrixData.at(i);
      }
      readRecord.geaneHits.geaneErrorMatrices.push_back(matrix);
    }
    readRecord.geaneHits.geaneErrorMatricesData.clear();

    for(auto & matrixData : record.geaneHits.extendedTransportMatrixBegToEndData){
      int matrixRows = sqrt(matrixData.size());
      int matrixCols = sqrt(matrixData.size());
      Eigen::MatrixXd matrix(matrixRows,matrixCols);
      for(unsigned int i = 0; i < matrixData.size(); i++){
        matrix(i) = matrixData.at(i);
      }
      readRecord.geaneHits.extendedTransportMatrixBegToEnd.push_back(matrix);
    }
    readRecord.geaneHits.extendedTransportMatrixBegToEndData.clear();


    // vector<Eigen::VectorXd>
    for(auto & vectorData : record.geaneHits.paramPredictedInUVEigenData){
      Eigen::VectorXd vector(vectorData.size());

      for(unsigned int i = 0; i < vectorData.size(); i++){
        vector(i) = vectorData.at(i);
      }
      readRecord.geaneHits.paramPredictedInUVEigen.push_back(vector);
    }
    readRecord.geaneHits.paramPredictedInUVEigenData.clear();


    return readRecord;
  }



  TrackArtRecord GeaneEigenStorageUtils::ReadEigenFromStorage(const TrackArtRecord& record)
  {
    // copy
    auto readRecord = record;

    // Eigen::MatrixXd
    int covarianceTotalInverseRows = sqrt(readRecord.geaneHits.covarianceTotalInverseData.size());
    int covarianceTotalInverseCols = sqrt(readRecord.geaneHits.covarianceTotalInverseData.size());
    readRecord.geaneHits.covarianceTotalInverse = Eigen::MatrixXd(covarianceTotalInverseRows,covarianceTotalInverseCols);
    for(unsigned int i = 0; i < readRecord.geaneHits.covarianceTotalInverseData.size(); i++){
       readRecord.geaneHits.covarianceTotalInverse(i) = readRecord.geaneHits.covarianceTotalInverseData.at(i);
    }
    readRecord.geaneHits.covarianceTotalInverseData.clear();

    for(auto & matrixData : record.geaneHits.geaneErrorMatricesData){
      int matrixRows = sqrt(matrixData.size());
      int matrixCols = sqrt(matrixData.size());
      Eigen::MatrixXd matrix(matrixRows,matrixCols);
      for(unsigned int i = 0; i < matrixData.size(); i++){
        matrix(i) = matrixData.at(i);
      }
      readRecord.geaneHits.geaneErrorMatrices.push_back(matrix);
    }
    readRecord.geaneHits.geaneErrorMatricesData.clear();

    for(auto & matrixData : record.geaneHits.extendedTransportMatrixBegToEndData){
      int matrixRows = sqrt(matrixData.size());
      int matrixCols = sqrt(matrixData.size());
      Eigen::MatrixXd matrix(matrixRows,matrixCols);
      for(unsigned int i = 0; i < matrixData.size(); i++){
        matrix(i) = matrixData.at(i);
      }
      readRecord.geaneHits.extendedTransportMatrixBegToEnd.push_back(matrix);
    }
    readRecord.geaneHits.extendedTransportMatrixBegToEndData.clear();

    return readRecord;
  }


  DecayVertexArtRecordCollection GeaneEigenStorageUtils::ReadEigenFromStorage(const DecayVertexArtRecordCollection& records)
  {
    // Declare collection that we're going to return (can't update in analyser as requires const argument)
    DecayVertexArtRecordCollection readCollection;

    // Loop over all records in collection
    for(auto & record : records){
      DecayVertexArtRecord readRecord = record;
      DecayVertexArtRecord updatedReadRecord = ReadEigenFromStorage(readRecord);          
      readCollection.push_back(updatedReadRecord);
    }

    // Return the updated collection
    return readCollection;
  }

  DecayVertexArtRecord GeaneEigenStorageUtils::ReadEigenFromStorage(const DecayVertexArtRecord& record)
  {
    // copy
    auto readRecord = record;

    // Eigen::MatrixXd
    int covarianceTotalRows = sqrt(readRecord.covMatrixData.size());
    int covarianceTotalCols = sqrt(readRecord.covMatrixData.size());
    readRecord.covMatrix = Eigen::MatrixXd(covarianceTotalRows,covarianceTotalCols);

    for(unsigned int i = 0; i < readRecord.covMatrixData.size(); i++){
       readRecord.covMatrix(i) = readRecord.covMatrixData.at(i);
    }
    readRecord.covMatrixData.clear();

    return readRecord;
  }

}//namespace
