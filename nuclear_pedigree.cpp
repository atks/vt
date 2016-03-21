#include "nuclear_pedigree.h"

NuclearPedigree::NuclearPedigree(const char* pedFile) {
  htsFile *hts = hts_open(pedFile, "r");
  kstring_t s = {0,0,0};
  std::vector<std::string> v;
  std::vector<std::string> smIDs;
  while(hts_getline(hts,'\n', &s)>=0) {
    if ( s.s[0] == '#' ) continue;
    std::string line(s.s);
    split(v, "\t", line);
    split(smIDs, ",", v[1]);
    addPerson(v[0].c_str(), smIDs, v[2] == "0" ? NULL : v[2].c_str(), v[3] == "0" ? NULL : v[3].c_str(), atoi(v[4].c_str()) );
  }
}

void NuclearPedigree::addPerson(const char* famID, const std::vector<std::string>& smIDs, const char* dadID, const char* momID, int sex) {

  // register family ID if absent
  NuclearFamily* pFam = NULL;
  if ( famIDmap.find(famID) == famIDmap.end() )  { // family ID does not exist
    pFam = new NuclearFamily(famID);
    famIDmap[famID] = pFam;
  }
  else 
    pFam = famIDmap[famID];

  // check whether the sample was seen before
  int i;
  for(i=0; i < (size_t)smIDs.size(); ++i) {
    if ( smIDmap.find(smIDs[i]) != smIDmap.end() )  { // sample ID has not seen before
      break;
    }
  }

  NuclearFamilyPerson* pPerson = NULL;
  if ( i < smIDs.size() ) {   // if sample was seen before, it must be dad or mom
    if ( ( pFam->pDad != NULL ) && ( pFam->pDad->hasSample(smIDs[i].c_str()) ) ) { // matches dad
      pPerson = pFam->pDad; // assign the person pointer
    }
    else if ( ( pFam->pMom != NULL ) && ( pFam->pMom->hasSample(smIDs[i].c_str()) ) ) {
      pPerson = pFam->pMom;
    }
    else {
      fprintf(stderr,"FATAL ERROR: Sample ID %s has already observed previously but the sample is not founders in the family %s\n", smIDs[i].c_str(), famID);
      exit(1);
    }
  }
  else {  // if the sample was not seen before, create the sample
    pPerson = new NuclearFamilyPerson(sex);
  }

  // register samples and assign to the person
  for(size_t j=0; j < smIDs.size(); ++j) {
    if ( smIDmap.find(smIDs[j]) == smIDmap.end() ) {  // register missing sample IDs
      NuclearFamilySample* p = new NuclearFamilySample(smIDs[j].c_str());
      smIDmap[smIDs[j]] = p;
      pPerson->addSample(p);
    }
    else { // if already registered, make sure that it belongs to the same person
      if ( !pPerson->hasSample(smIDs[j].c_str()) ) {
	fprintf(stderr,"FATAL ERROR: Sample ID %s has already observed previously but the sample is not founders in the family %s\n", smIDs[j].c_str(), famID);
	exit(1);	
      }
    }
  }

  // add the person to the family
  if ( ( dadID == NULL ) && ( momID == NULL ) ) { // founders
    if ( sex == 1 ) {  // dad
      if ( ( pFam->pDad != NULL ) && ( pFam->pDad != pPerson ) ) {
	fprintf(stderr,"FATAL ERROR: Multiple fathers %s, %s are observed for family %s\n", pPerson->getUID(), pFam->pDad->getUID(), famID);
	exit(1);
      }
      pFam->pDad = pPerson;
    }
    else if ( sex == 2 ) {
      if ( ( pFam->pMom != NULL ) && ( pFam->pMom != pPerson ) ) {
	fprintf(stderr,"FATAL ERROR: Multiple fathers %s, %s are observed for family %s\n", pPerson->getUID(), pFam->pMom->getUID(), famID);
	exit(1);
      }
      pFam->pMom = pPerson;
    }
    else {
      fprintf(stderr,"FATAL ERROR: Founder %d must have sex information available", smIDs[0].c_str());
      exit(1);
    }
  }
  else { // not a founder - then add as children
    pFam->pKids.push_back(pPerson);
  }
}

bool NuclearPedigree::setSampleIndex(const char* smID, int _index) {
  if ( smIDmap.find(smID) == smIDmap.end() )
    return false;
  else {
    smIDmap[smID]->index = _index;
    return true;
  }
}

int NuclearPedigree::numSamplesWithIndex() {
  int count = 0;
  std::map<std::string, NuclearFamilySample*>::iterator it;
  for(it = smIDmap.begin(); it != smIDmap.end(); ++it) {
    if ( it->second->index >= 0 )
      ++count;
  }
  return count;
}

int NuclearPedigree::numPeople() {
  int count = 0;
  std::map<std::string, NuclearFamily*>::iterator it;
  for(it = famIDmap.begin(); it != famIDmap.end(); ++it) {
    NuclearFamily* pFam = it->second;    
    if ( pFam->pDad ) ++count;
    if ( pFam->pMom ) ++count;
    count += (int)pFam->pKids.size();
  }
  return count;
}

int NuclearPedigree::removeSamplesWithoutIndex() {
  std::map<std::string, NuclearFamily*>::iterator it;
  int nRemoved = 0;
  for(it = famIDmap.begin(); it != famIDmap.end(); ++it) {
    NuclearFamily* pFam = it->second;
    if ( pFam->pDad ) {
      int org = (int)pFam->pDad->samples.size();
      int rm = pFam->pDad->removeSamplesWithoutIndex(smIDmap);
      if ( org == rm ) {
	delete pFam->pDad;
	pFam->pDad = NULL;
      }
      nRemoved += rm;
    }
    if ( pFam->pMom ) {
      int org = (int)pFam->pMom->samples.size();
      int rm = pFam->pMom->removeSamplesWithoutIndex(smIDmap);
      if ( org == rm ) {      
	delete pFam->pMom;
	pFam->pMom = NULL;
      }
      nRemoved += rm;      
    }
    for(int i=(int)pFam->pKids.size()-1; i>=0; --i) {
      int org = (int)pFam->pKids[i]->samples.size();
      int rm = pFam->pKids[i]->removeSamplesWithoutIndex(smIDmap);
      if ( org == rm ) {            
	delete pFam->pKids[i];
	pFam->pKids.erase(pFam->pKids.begin()+i);
      }
      nRemoved += rm;            
    }

    // remove family if no sample exists in the VCF/BCF file
    if (pFam->pKids.empty() && pFam->pDad == NULL && pFam->pMom == NULL) {
      famIDmap.erase(pFam->famID);
      delete pFam;
    }
  }
  return nRemoved;
}
