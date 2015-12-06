/**
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include "btree.h"
#include "filescan.h"
#include "exceptions/bad_index_info_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/end_of_file_exception.h"
#include "exceptions/file_exists_exception.h"
#include "exceptions/hash_not_found_exception.h"
#include "exceptions/page_not_pinned_exception.h"


//#define DEBUG

namespace badgerdb
{

// -----------------------------------------------------------------------------
// BTreeIndex::BTreeIndex -- Constructor
// -----------------------------------------------------------------------------

BTreeIndex::BTreeIndex(const std::string & relationName,
		std::string & outIndexName,
		BufMgr *bufMgrIn,
		const int attrByteOffset,
		const Datatype attrType)
{
	std::ostringstream idxStr;
	idxStr << relationName << '.' << attrByteOffset;
	outIndexName = idxStr.str();

	this->bufMgr = bufMgrIn;
	this->attrByteOffset = attrByteOffset;
	this->attributeType = attrType;
	scanExecuting = false;
	leafOccupancy = 0;
	nodeOccupancy = 0;

	try{
		//metadata page
		file = new BlobFile(outIndexName, true);
		Page* page;
		bufMgrIn->allocPage(file, headerPageNum, page);
		IndexMetaInfo* indexMetaInfo = (IndexMetaInfo*) page;
		indexMetaInfo->attrByteOffset = attrByteOffset;
		indexMetaInfo->attrType = attrType;
		strcpy(indexMetaInfo->relationName, relationName.c_str());
		//root page
		bufMgrIn->allocPage(file, rootPageNum, page);
		indexMetaInfo->rootPageNo = rootPageNum;
		bufMgrIn->unPinPage(file, headerPageNum, true);

		if(attributeType == INTEGER){
			NonLeafNodeInt *nonLeafNodeInt = (NonLeafNodeInt *) page;
			for(int i=0;i<INTARRAYNONLEAFSIZE;i++) nonLeafNodeInt->pageNoArray[i] = 0;
			nonLeafNodeInt->level = 1;
		} else if(attributeType == DOUBLE){
			NonLeafNodeDouble *nonLeafNodeDouble = (NonLeafNodeDouble *) page;
			for(int i=0;i<DOUBLEARRAYNONLEAFSIZE;i++) nonLeafNodeDouble->pageNoArray[i] = 0;
			nonLeafNodeDouble->level = 1;
		} else {
			NonLeafNodeString *nonLeafNodeString = (NonLeafNodeString *) page;
			for(int i=0;i<STRINGARRAYNONLEAFSIZE;i++) nonLeafNodeString->pageNoArray[i] = 0;
			nonLeafNodeString->level = 1;
		}
		bufMgrIn->unPinPage(file, rootPageNum, true);

		//insert index
		FileScan fscan(relationName, bufMgr);
		RecordId scanRid;
		while(1)
		{
			fscan.scanNext(scanRid);
			insertEntry((void *) fscan.getRecord().c_str() + attrByteOffset, scanRid);
		}
	}
	catch (FileExistsException e){
		this->file = new BlobFile(outIndexName, false);
		Page* page;
		headerPageNum = 1;
		bufMgrIn->readPage(file, headerPageNum, page);
		IndexMetaInfo* indexMetaInfo = (IndexMetaInfo*) page;

		//check
		if(indexMetaInfo->attrByteOffset != attrByteOffset
		   || indexMetaInfo->attrType != attrType
		   || strcmp(indexMetaInfo->relationName, relationName.c_str()) != 0
				){
			bufMgrIn->unPinPage(file, headerPageNum, false);
			throw BadIndexInfoException("constructor parameters do not match exist index file");
		}

		this->rootPageNum = indexMetaInfo->rootPageNo;
		bufMgrIn->unPinPage(file, headerPageNum, false);

		std::cout << "Finish Read index file "
		<< " hearpage: " << headerPageNum << " rootpage: " << rootPageNum<< std::endl; //TODO: delete
	}
	catch (EndOfFileException e){
		std::cout << "Finish Read all records "
		<< " hearpage: " << headerPageNum << " rootpage: " << rootPageNum<< std::endl; //TODO: delete
	}
}


// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------

BTreeIndex::~BTreeIndex()
{
	scanExecuting = false;

	try {
		bufMgr->unPinPage(file, currentPageNum, false);
	} catch (PageNotPinnedException e) {
	} catch (HashNotFoundException e) {
	}

	bufMgr->flushFile(file);
	file->~File();
}

template <class T>
void copy(T& a, T& b){
		a = b;
}

template <>
void copy<char[STRINGSIZE]>(char (&a)[STRINGSIZE], char (&b)[STRINGSIZE]){
	strncpy( a, b, STRINGSIZE);
	//a[STRINGSIZE - 1] = '\0';
}

template <class T>
bool smallerThan(T a, T b){
	return a<b;
}

template <>
bool smallerThan<char[STRINGSIZE]>(char a[STRINGSIZE], char b[STRINGSIZE]){
	return strncmp(a, b, STRINGSIZE) < 0;
}

template <class T>
bool biggerThan(T a, T b){
	return a>b;
}

template <>
bool biggerThan<char[STRINGSIZE]>(char a[STRINGSIZE], char b[STRINGSIZE]){
	return strncmp(a, b, STRINGSIZE) > 0;		
}

template <class T>
bool equalTo(T a, T b){
	return a==b;
}

template <>
bool equalTo<char[STRINGSIZE]>(char a[STRINGSIZE], char b[STRINGSIZE]){
	return strncmp(a, b, STRINGSIZE) == 0;
}



template <class T, class T1>
void BTreeIndex::leafSplitHelper(int pos, int last, int LEAFARRAYMAX,
								 int NONLEAFARRAYMAX,
								 RIDKeyPair<T> ridKeyPair,
								 T1* leafNode,
							  	 PageId& newPageId,
								 T & newValue)
{
	//full
	Page* newPage;
	bufMgr->allocPage(file, newPageId, newPage);
	T1 *newLeafNode = (T1 *) newPage;

	//tmp array
	T tmpKeyArray[LEAFARRAYMAX + 1];
	RecordId tmpRidArray[LEAFARRAYMAX + 1];

	//copy all new records to tmp
	for(int i=0; i < LEAFARRAYMAX; i++) {
		tmpRidArray[i] = leafNode->ridArray[i];
		copy<T> (tmpKeyArray[i], leafNode->keyArray[i]);
		leafNode->ridArray[i].page_number = 0;
	}

	for(int i= LEAFARRAYMAX; i > pos; i--){
		copy<T> (tmpKeyArray[i], tmpKeyArray[i-1]);
		tmpRidArray[i] = tmpRidArray[i-1];
	}
	copy<T> (tmpKeyArray[pos], ridKeyPair.key);
	tmpRidArray[pos] = ridKeyPair.rid;

	//reset new and old page
	for(int i=0; i < LEAFARRAYMAX; i++) {
		leafNode->ridArray[i].page_number = 0;
		newLeafNode->ridArray[i].page_number = 0;
	}

	//copy back
	for(int i=0; i< (LEAFARRAYMAX + 1) / 2; i++ ){
		copy<T> (leafNode->keyArray[i], tmpKeyArray[i]);
		leafNode->ridArray[i] = tmpRidArray[i];
	}
	for(int i= (LEAFARRAYMAX + 1) / 2; i < LEAFARRAYMAX + 1; i++){
		copy<T> (newLeafNode->keyArray[i - (LEAFARRAYMAX + 1) / 2], tmpKeyArray[i]);
		newLeafNode->ridArray[i - (LEAFARRAYMAX + 1) / 2] = tmpRidArray[i];
	}

	//link leaf node
	newLeafNode->rightSibPageNo = leafNode->rightSibPageNo;
	leafNode->rightSibPageNo = newPageId;

	//push up
	copy<T> (newValue, newLeafNode->keyArray[0]);

	//unpin
	bufMgr->unPinPage(file, newPageId, true);
};



// insertEntry when is full
template <class T, class T2>
void  BTreeIndex::nonLeafSplitHelper(int pos,
									 int NONLEAFARRAYMAX,
									 T2* nonLeafNode,
									 PageId& newPageId,
									 T & newValue,
									 T&  newChildValue,
									 PageId newChildPageId){
	//full, need split
	Page* newPage;
	bufMgr->allocPage(file, newPageId, newPage);
	T2* newNonLeafNode = (T2*) newPage;

	//tmp array
	T tmpKeyArray[NONLEAFARRAYMAX + 1];
	PageId tmpPageIdArray[NONLEAFARRAYMAX + 2];

	//copy to tmp
	for(int i=0; i < NONLEAFARRAYMAX; i++) {
		tmpPageIdArray[i] = nonLeafNode->pageNoArray[i];
		copy<T> (tmpKeyArray[i], nonLeafNode->keyArray[i]);
	}
	tmpPageIdArray[NONLEAFARRAYMAX + 1] = nonLeafNode->pageNoArray[NONLEAFARRAYMAX + 1];

	for(int i= NONLEAFARRAYMAX; i > pos; i--){
		copy<T> (tmpKeyArray[i], tmpKeyArray[i-1]);
		tmpPageIdArray[i+1] = tmpPageIdArray[i];
	}
	copy<T> (tmpKeyArray[pos], newChildValue);
	tmpPageIdArray[pos+1] = newChildPageId;

	//clear old and new page
	for(int i=0; i < NONLEAFARRAYMAX + 1; i++){
		nonLeafNode->pageNoArray[i] = 0;
		newNonLeafNode->pageNoArray[i] = 0;
	}

	//copy back
	for(int i=0; i< (NONLEAFARRAYMAX + 1) / 2; i++ ){
		copy<T> (nonLeafNode->keyArray[i], tmpKeyArray[i]);
		nonLeafNode->pageNoArray[i] = tmpPageIdArray[i];
	}
	nonLeafNode->pageNoArray[(NONLEAFARRAYMAX + 1) / 2] = tmpPageIdArray[(NONLEAFARRAYMAX + 1) / 2];

	for(int i= (NONLEAFARRAYMAX + 1) / 2 + 1; i < NONLEAFARRAYMAX + 1; i++){
		copy<T> (newNonLeafNode->keyArray[i - (NONLEAFARRAYMAX + 1) / 2 - 1], tmpKeyArray[i]);
		newNonLeafNode->pageNoArray[i - (NONLEAFARRAYMAX + 1) / 2 - 1 ] = tmpPageIdArray[i];
	}
	newNonLeafNode->pageNoArray[NONLEAFARRAYMAX + 1 - (NONLEAFARRAYMAX + 1) / 2 - 1 ] = tmpPageIdArray[
			NONLEAFARRAYMAX + 1];

	//level
	newNonLeafNode->level = nonLeafNode->level;

	//push up
	copy<T> (newValue, tmpKeyArray[(NONLEAFARRAYMAX + 1) / 2]);

	//unpin
	bufMgr->unPinPage(file, newPageId, true);
};



template <class T, class T1, class T2>
void BTreeIndex::createFirstLeaf(int LEAFARRAYMAX,
								 RIDKeyPair<T> ridKeyPair,
								 T2* nonLeafNode,
								 PageId pageId)
{
	std::cout<<"empty index"  << std::endl; //TODO:delete
	PageId newPageId;
	Page* page;
	bufMgr->allocPage(file, newPageId, page);

	T1* leafNode = (T1*) page;
	for(int i=0; i < LEAFARRAYMAX; i++)
		leafNode->ridArray[i].page_number = 0;

	copy<T> (leafNode->keyArray[0], ridKeyPair.key);
	leafNode->ridArray[0] = ridKeyPair.rid;
	leafNode->rightSibPageNo = 0;
	nonLeafNode->pageNoArray[0] = newPageId;

	//unpin
	bufMgr->unPinPage(file, newPageId, true);
	bufMgr->unPinPage(file, pageId, true);

};





template <class T, class T1, class T2>
void BTreeIndex::insertEntryRecursive(RIDKeyPair<T> ridKeyPair,
									  PageId pageId,
									  bool isLeaf,
									  int LEAFARRAYMAX,
									  int NONLEAFARRAYMAX,
									  T & newValue,
									  PageId& newPageId)
{
//	std::cout<<"try to insert key: "<< ridKeyPair.key << std::endl; //TODO:delete

	if(isLeaf){
		//read the page
		Page* page;
		bufMgr->readPage(file, pageId, page);
		T1 *leafNode = (T1 *) page;

		//find position
		int pos = 0;
		while(biggerThan<T> (ridKeyPair.key, leafNode->keyArray[pos])
			  && leafNode->ridArray[pos].page_number != 0
			  && pos < LEAFARRAYMAX)
			pos++;

		//find last entry
		int last;
		for (last =0; last < LEAFARRAYMAX; last++)
			if(leafNode->ridArray[last].page_number == 0)
				break;


		if(last < LEAFARRAYMAX){
			//not full
			for(int i=last; i>pos; i--){
				copy<T> (leafNode->keyArray[i], leafNode->keyArray[i - 1]);
				leafNode->ridArray[i] = leafNode->ridArray[i - 1];
			}
			copy<T> (leafNode->keyArray[pos], ridKeyPair.key);
			leafNode->ridArray[pos] = ridKeyPair.rid;

			std::cout<<"insert key: "<< ridKeyPair.key <<" at page :" << pageId
					<< " at position: "<< pos <<std::endl; //TODO:delete
		} else {
			leafSplitHelper<T,T1>(pos, last, LEAFARRAYMAX,
					NONLEAFARRAYMAX, ridKeyPair,leafNode,newPageId,newValue);
//			std::cout<<"insert key: "<< ridKeyPair.key <<" at page :" << newPageId << std::endl; //TODO:delete
		}

		//unpin
		bufMgr->unPinPage(file, pageId, true);

	} else {
		//non leaf
		//read page
		Page* page;
		bufMgr->readPage(file, pageId, page);
		T2* nonLeafNode = (T2*) page;

		//find pageArray position
		int pos = 0;
		
		while(!smallerThan<T> (ridKeyPair.key, nonLeafNode->keyArray[pos])
			  && nonLeafNode->pageNoArray[pos + 1] != 0
			  && pos < NONLEAFARRAYMAX)
			pos++;

		//index file is empty
		if(nonLeafNode->pageNoArray[pos] == 0){
			createFirstLeaf<T, T1, T2>(LEAFARRAYMAX,
					ridKeyPair,
					nonLeafNode,
					pageId);
			return;
		}

		//call recursive function
		bufMgr->unPinPage(file, pageId, false);
		T newChildValue;
		PageId newChildPageId = 0;
		insertEntryRecursive<T, T1, T2>(ridKeyPair, nonLeafNode->pageNoArray[pos], nonLeafNode->level == 1,
										LEAFARRAYMAX, NONLEAFARRAYMAX, newChildValue, newChildPageId);

		//check if child split
		if(newChildPageId != 0){
			//if child split
			bufMgr->readPage(file, pageId, page);
			T2* nonLeafNode = (T2*) page;

			//find last entry
			int last;
			for (last = 0; last < NONLEAFARRAYMAX; last++)
				if(nonLeafNode->pageNoArray[last + 1] == 0)
					break;

			//check if full
			if(last < NONLEAFARRAYMAX){
				//not full
				for(int i=last; i>pos; i--){
					copy<T> (nonLeafNode->keyArray[i], nonLeafNode->keyArray[i - 1]);
					nonLeafNode->pageNoArray[i + 1] = nonLeafNode->pageNoArray[i];
				}
				copy<T> (nonLeafNode->keyArray[pos], newChildValue);
				nonLeafNode->pageNoArray[pos + 1] = newChildPageId;
			}
			else{
				nonLeafSplitHelper<T, T2>(pos, NONLEAFARRAYMAX, nonLeafNode,newPageId,newValue,
										  newChildValue,newChildPageId);
			}

			bufMgr->unPinPage(file, pageId, true);
		}

		//unpin
		//bufMgr->unPinPage(file, pageId, true);
	}
}

template<class T, class T1>
void BTreeIndex::handleNewRoot(T& newValue, PageId newPageId, int ARRAYMAX){
	std::cout<<"handle new root"<<std::endl; //TODO
	PageId newRootPageId;
	Page *newRootPage;
	bufMgr->allocPage(file, newRootPageId, newRootPage);

	T1 *newRootNonLeafNode = (T1 *) newRootPage;
	for(int i=0;i<ARRAYMAX + 1; i++) newRootNonLeafNode->pageNoArray[i] = 0;
	copy<T> (newRootNonLeafNode->keyArray[0], newValue);
	newRootNonLeafNode->pageNoArray[0] = rootPageNum;
	newRootNonLeafNode->pageNoArray[1] = newPageId;
	newRootNonLeafNode->level = 0;
	rootPageNum = newRootPageId;
	bufMgr->unPinPage(file, newRootPageId, true);std::cout<<"handle new root finish"<<std::endl; //TODO
	std::cout<<"handle new root finish"<<std::endl; //TODO
}
// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

const void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{

	if(attributeType == INTEGER) {

		RIDKeyPair<int> ridKeyPairInt;
		ridKeyPairInt.set(rid, *((int *) key));
		int newValueInt;
		PageId newPageIdInt = 0;

		//call recursive function
		insertEntryRecursive<int, LeafNodeInt, NonLeafNodeInt>
				(ridKeyPairInt, rootPageNum, 0, INTARRAYLEAFSIZE, INTARRAYNONLEAFSIZE, newValueInt, newPageIdInt);

		//if root got split
		if (newPageIdInt != 0)
			handleNewRoot<int, NonLeafNodeInt>(newValueInt, newPageIdInt, INTARRAYNONLEAFSIZE);
	} else if (attributeType == DOUBLE) {

		RIDKeyPair<double> ridKeyPairDouble;
		ridKeyPairDouble.set(rid, *((double *) key));
		double newValueDouble;
		PageId newPageIdDouble = 0;

		//call recursive function
		insertEntryRecursive<double, LeafNodeDouble, NonLeafNodeDouble>
				(ridKeyPairDouble, rootPageNum, 0, DOUBLEARRAYLEAFSIZE, DOUBLEARRAYNONLEAFSIZE, newValueDouble, newPageIdDouble);

		//if root got split
		if (newPageIdDouble != 0)
			handleNewRoot<double, NonLeafNodeDouble>(newValueDouble, newPageIdDouble, DOUBLEARRAYNONLEAFSIZE);

	} else if (attributeType == STRING){

		RIDKeyPair<char[STRINGSIZE] > ridKeyPairString;
		ridKeyPairString.rid = rid;
		strncpy(ridKeyPairString.key, (char*) key, STRINGSIZE);
		//ridKeyPairString.key[STRINGSIZE - 1] = '\0';
		char newValue[STRINGSIZE];
		PageId newPageId = 0;

		//call recursive function
		insertEntryRecursive<char[STRINGSIZE], LeafNodeString, NonLeafNodeString>
				(ridKeyPairString, rootPageNum, 0, STRINGARRAYLEAFSIZE, STRINGARRAYNONLEAFSIZE, newValue, newPageId);

		//if root got split
		if (newPageId != 0)
			handleNewRoot<char[STRINGSIZE], NonLeafNodeString>(newValue, newPageId, STRINGARRAYNONLEAFSIZE);

	}

}

// -----------------------------------------------------------------------------
// BTreeIndex::startScan
// -----------------------------------------------------------------------------

const void BTreeIndex::startScan(const void* lowValParm,
								 const Operator lowOpParm,
								 const void* highValParm,
								 const Operator highOpParm)
{
	scanExecuting = true;
	if (lowOpParm != GT && lowOpParm !=GTE ){
		throw BadOpcodesException();
	}
	if(highOpParm != LT && highOpParm != LTE){
		throw BadOpcodesException();
	}

	this->lowOp = lowOpParm;
	this->highOp = highOpParm;

	if(attributeType == INTEGER) {
		this->lowValInt = *((int*) lowValParm);
		this->highValInt = *((int *) highValParm);
		startScanHelper<int, NonLeafNodeInt>(*((int*) lowValParm), *((int *) highValParm), INTARRAYNONLEAFSIZE);
	} else if (attributeType == DOUBLE) {
		this->lowValDouble = *((double*) lowValParm);
		this->highValDouble = *((double *) highValParm);
		startScanHelper<double , NonLeafNodeDouble>(*((double *) lowValParm), *((double *) highValParm), DOUBLEARRAYNONLEAFSIZE);
	} else {
		strncpy((char*) this->lowValString.c_str(), (char *)lowValParm, STRINGSIZE-1);
		strncpy(lowValChar, (char *)lowValParm, STRINGSIZE);
		//lowValChar[STRINGSIZE - 1] = '\0';
		strncpy((char*) this->highValString.c_str(), (char *)highValParm, STRINGSIZE-1);
		strncpy(highValChar, (char *)highValParm, STRINGSIZE);
		//highValChar[STRINGSIZE - 1] = '\0';

		std::cout<< "Start scan finish, lowValChar: "<< lowValChar <<" highValChar: "<< highValChar<<std::endl;
		startScanHelper<char[STRINGSIZE], NonLeafNodeString> (lowValChar, highValChar, STRINGARRAYNONLEAFSIZE);
	}

}

template<class T, class T1>
void BTreeIndex::startScanHelper(T lowValParm,
								 T highValParm,
								 int NONLEAFARRAYMAX)
{
	if(biggerThan<T> (lowValParm,  highValParm))
		throw BadScanrangeException();

	//find first one
	currentPageNum = rootPageNum;	std::cout<<"start scan read root page:" << currentPageNum<<std::endl;//TODO
	bufMgr->readPage(file, currentPageNum, currentPageData);
	T1* nonLeafNode = (T1*) currentPageData; std::cout<<"current page level:" << nonLeafNode->level<<std::endl;//TODO

	int pos = 0;
	while(nonLeafNode->level != 1) {std::cout<<"start scan go to next nonleaf page:" << nonLeafNode->pageNoArray[0]<<std::endl;//TODO
		pos = 0;
		while(!smallerThan<T> (lowValParm, nonLeafNode->keyArray[pos])
			  && nonLeafNode->pageNoArray[pos + 1] != 0
			  && pos < NONLEAFARRAYMAX)
			pos++;
		PageId nextPageId = nonLeafNode->pageNoArray[pos];
		bufMgr->readPage(file, nextPageId, currentPageData);
		bufMgr->unPinPage(file, currentPageNum, false);
		currentPageNum = nextPageId;
		nonLeafNode = (T1*) currentPageData;
	}
	pos = 0;
	while(!smallerThan<T> (lowValParm, nonLeafNode->keyArray[pos])
		  && nonLeafNode->pageNoArray[pos + 1] != 0
		  && pos < NONLEAFARRAYMAX)
		pos++;
	PageId nextPageId = nonLeafNode->pageNoArray[pos];
	bufMgr->readPage(file, nextPageId, currentPageData);
	bufMgr->unPinPage(file, currentPageNum, false);
	currentPageNum = nextPageId;
	std::cout<<"start scan go to next leaf page:" << currentPageNum <<std::endl;//TODO
	nextEntry = 0;

	std::cout<< "start scan complete" << std::endl;
}


// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

const void BTreeIndex::scanNext(RecordId& outRid) 
{
	if(scanExecuting == false)
		throw ScanNotInitializedException();

	if(attributeType == INTEGER)
		scanNextHelper<int, LeafNodeInt> (outRid, lowValInt, highValInt, INTARRAYLEAFSIZE);
	else if (attributeType == DOUBLE)
		scanNextHelper<double, LeafNodeDouble> (outRid, lowValDouble, highValDouble, DOUBLEARRAYLEAFSIZE);
	else
		scanNextHelper<char[STRINGSIZE], LeafNodeString> (outRid, lowValChar, highValChar, STRINGARRAYLEAFSIZE);
}

template <class T, class T1>
void BTreeIndex::scanNextHelper(RecordId &outRid, T lowVal, T highVal, int ARRAYMAX)
{
	T1* leafNode;
	while(1){
		leafNode = (T1*) currentPageData;

		//Go to next page
		if(leafNode->ridArray[nextEntry].page_number == 0 || nextEntry == ARRAYMAX) {

			PageId nextPageNum = leafNode->rightSibPageNo;
			if(nextPageNum == 0){
				std::cout<<"scan finish"<<std::endl; //TODO:delete
				bufMgr->unPinPage(file, currentPageNum, false);
				throw IndexScanCompletedException();
			}

			bufMgr->unPinPage(file, currentPageNum, false);
			currentPageNum = nextPageNum;

			bufMgr->readPage(file, currentPageNum, currentPageData);
			std::cout<<std::endl <<"scan next go to next leaf page: "<< nextPageNum <<std::endl; //TODO:delete
			nextEntry = 0;
			continue;
		}
//		std::cout<<" Got: "<< leafNode->keyArray[nextEntry] ; //TODO:delete
		//Do not satisfy
		if((lowOp==GT && !biggerThan<T> (leafNode->keyArray[nextEntry], lowVal) )
		   || (lowOp==GTE && smallerThan<T> (leafNode->keyArray[nextEntry], lowVal))
		   ) {
			nextEntry++;
			continue;
		}

		if((highOp==LT && !smallerThan<T> (leafNode->keyArray[nextEntry],  highVal))
		   || (highOp==LTE && biggerThan<T> (leafNode->keyArray[nextEntry], highVal) ))
			throw IndexScanCompletedException();


//		std::cout<<" Got: "<< leafNode->keyArray[nextEntry] ; //TODO:delete


		//Got a record
		outRid = leafNode->ridArray[nextEntry];
		nextEntry++;
		return ;
	}
}
// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//
const void BTreeIndex::endScan() 
{
	if(scanExecuting == false)
		throw ScanNotInitializedException();
	scanExecuting = false;

	try {
		bufMgr->unPinPage(file, currentPageNum, false);
	} catch (PageNotPinnedException e) {
	} catch (HashNotFoundException e) {
	}

//	bufMgr->flushFile(file);
}


}
