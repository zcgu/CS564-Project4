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

	try{
		//metadata page
		file = new BlobFile(outIndexName, true);
		Page* page;
		bufMgrIn->allocPage(file, headerPageNum, page);
		IndexMetaInfo* indexMetaInfo = (IndexMetaInfo*) page;
		indexMetaInfo->attrByteOffset = attrByteOffset;
		indexMetaInfo->attrType = attrType;
		strcpy(indexMetaInfo->relationName, relationName.c_str());
		bufMgrIn->unPinPage(file, headerPageNum, true);

		//root page
		bufMgrIn->allocPage(file, rootPageNum, page);
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

		bufMgrIn->readPage(file, headerPageNum, page);
		indexMetaInfo = (IndexMetaInfo*) page;
		indexMetaInfo->rootPageNo = rootPageNum;
		bufMgrIn->unPinPage(file, headerPageNum, true);

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
				)
			throw BadIndexInfoException("constructor parameters do not match exist index file");

		this->rootPageNum = indexMetaInfo->rootPageNo;
		bufMgrIn->unPinPage(file, headerPageNum, false);

		std::cout << "Finish Read index file //Gu"
		<< "hearpage: " << headerPageNum << "rootpage:" << rootPageNum<< std::endl; //TODO: delete
	}
	catch (EndOfFileException e){
		std::cout << "Finish Read all records //Gu"
		<< "hearpage: " << headerPageNum << "rootpage:" << rootPageNum<< std::endl; //TODO: delete
	}
}


// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------

BTreeIndex::~BTreeIndex()
{
	bufMgr->flushFile(file);
	file->~File();
}


/*
void BTreeIndex::insertEntryRecursive(RIDKeyPair<T > ridKeyPair,
									  PageId pageId,
									  bool isLeaf,
						  				T& newValue,
						  			PageId& newPageId)
{
	if(isLeaf){
		Page* page;
		bufMgr->readPage(file, pageId, page);
		switch (attributeType){
			case INTEGER:{
				LeafNodeInt* leafNodeInt = (LeafNodeInt*) page;
				int pos = 0;
				RIDKeyPair<T> tmpRidKeyPair;
				RecordId tmpRecordId;
				tmpRidKeyPair.set(tmpRecordId, leafNodeInt->keyArray[pos]);
				while( !(ridKeyPair < tmpRidKeyPair) ){
					if( ++pos == INTARRAYLEAFSIZE) break;
					tmpRidKeyPair.set(tmpRecordId, leafNodeInt->keyArray[pos]);
				}

				int last;
				for (last =0; last < INTARRAYLEAFSIZE; last++)
					if(leafNodeInt->ridArray[last].page_number == 0)
						break;

				if(last < INTARRAYLEAFSIZE){
					//not full
					for(int i=last; i>pos; i--){
						leafNodeInt->keyArray[i] = leafNodeInt->keyArray[i-1];
						leafNodeInt->ridArray[i] = leafNodeInt->ridArray[i-1];
					}
					leafNodeInt->keyArray[pos] = ridKeyPair.key;
					leafNodeInt->ridArray[pos] = ridKeyPair.rid;
				} else {
					//full
					Page* newPage;
					bufMgr->allocPage(file, newPageId, newPage);

					//tmp array
					int tmpKeyArray[INTARRAYLEAFSIZE+1];
					RecordId tmpRidArray[INTARRAYLEAFSIZE +1];

					//copy to tmp
					for(int i=0; i<INTARRAYLEAFSIZE; i++) {
						tmpRidArray[i] = leafNodeInt->ridArray[i];
						tmpKeyArray[i] = leafNodeInt->keyArray[i];
						leafNodeInt->ridArray[i].page_number = 0;
					}

					for(int i=INTARRAYLEAFSIZE; i>pos; i--){
						tmpKeyArray[i] = tmpKeyArray[i-1];
						tmpRidArray[i] = tmpRidArray[i-1];
					}
					tmpKeyArray[pos] = ridKeyPair.key;
					tmpRidArray[pos] = ridKeyPair.rid;

					//copy back
					for(int i=0; i<(INTARRAYLEAFSIZE+1)/2;i++ ){
;
					}


				}

			}
			case DOUBLE:{

			}
			case STRING:{

			}
		}
	} else {
		//non leaf


	}
}
*/

template <class T, class T1, class T2>
void BTreeIndex::insertEntryRecursive(RIDKeyPair<T> ridKeyPair,
									  PageId pageId,
									  bool isLeaf,
									  int ARRAYMAX1,
									  int ARRAYMAX2,
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
		while(ridKeyPair.key > leafNode->keyArray[pos]
			  && leafNode->ridArray[pos].page_number != 0
			  && pos < ARRAYMAX1)
			pos++;

		//find last entry
		int last;
		for (last =0; last < ARRAYMAX1; last++)
			if(leafNode->ridArray[last].page_number == 0)
				break;


		if(last < ARRAYMAX1){
			//not full
			for(int i=last; i>pos; i--){
				leafNode->keyArray[i] = leafNode->keyArray[i - 1];
				leafNode->ridArray[i] = leafNode->ridArray[i - 1];
			}
			leafNode->keyArray[pos] = ridKeyPair.key;
			leafNode->ridArray[pos] = ridKeyPair.rid;

//			std::cout<<"insert key: "<< ridKeyPair.key <<" at page :" << pageId << std::endl; //TODO:delete
		} else {
			//full
			Page* newPage;
			bufMgr->allocPage(file, newPageId, newPage);
			T1 *newLeafNode = (T1 *) newPage;

			//tmp array
			int tmpKeyArray[ARRAYMAX1 + 1];
			RecordId tmpRidArray[ARRAYMAX1 + 1];

			//copy all new records to tmp
			for(int i=0; i < ARRAYMAX1; i++) {
				tmpRidArray[i] = leafNode->ridArray[i];
				tmpKeyArray[i] = leafNode->keyArray[i];
				leafNode->ridArray[i].page_number = 0;
			}

			for(int i= ARRAYMAX1; i > pos; i--){
				tmpKeyArray[i] = tmpKeyArray[i-1];
				tmpRidArray[i] = tmpRidArray[i-1];
			}
			tmpKeyArray[pos] = ridKeyPair.key;
			tmpRidArray[pos] = ridKeyPair.rid;

			//reset new and old page
			for(int i=0; i < ARRAYMAX1; i++) {
				leafNode->ridArray[i].page_number = 0;
				newLeafNode->ridArray[i].page_number = 0;
			}

			//copy back
			for(int i=0; i< (ARRAYMAX1 + 1) / 2; i++ ){
				leafNode->keyArray[i] = tmpKeyArray[i];
				leafNode->ridArray[i] = tmpRidArray[i];
			}
			for(int i= (ARRAYMAX1 + 1) / 2; i < ARRAYMAX1 + 1; i++){
				newLeafNode->keyArray[i - (ARRAYMAX1 + 1) / 2] = tmpKeyArray[i];
				newLeafNode->ridArray[i - (ARRAYMAX1 + 1) / 2] = tmpRidArray[i];
			}

			//link leaf node
			newLeafNode->rightSibPageNo = leafNode->rightSibPageNo;
			leafNode->rightSibPageNo = newPageId;

			//push up
			newValue = newLeafNode->keyArray[0];

			//unpin
			bufMgr->unPinPage(file, newPageId, true);

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
		while(ridKeyPair.key >= nonLeafNode->keyArray[pos]
			  && nonLeafNode->pageNoArray[pos + 1] != 0
			  && pos < ARRAYMAX2)
			pos++;

		//index file is empty
		if(nonLeafNode->pageNoArray[pos] == 0){
			std::cout<<"empty index"  << std::endl; //TODO:delete
			PageId newPageIdLeft, newPageIdRight;
			Page* pageLeft, *pageRight;
			bufMgr->allocPage(file, newPageIdLeft, pageLeft);
			bufMgr->allocPage(file, newPageIdRight, pageRight);

			T1* leafNodeLeft = (T1*) pageLeft;
			T1* leafNodeRight = (T1*) pageRight;
			for(int i=0; i < ARRAYMAX2; i++){
				leafNodeLeft->ridArray[i].page_number = 0;
				leafNodeRight->ridArray[i].page_number = 0;
			}
			leafNodeRight->keyArray[0] = ridKeyPair.key;
			leafNodeRight->ridArray[0] = ridKeyPair.rid;
			leafNodeLeft->rightSibPageNo = newPageIdRight;
			leafNodeRight->rightSibPageNo = 0;
			nonLeafNode->keyArray[0] = ridKeyPair.key;
			nonLeafNode->pageNoArray[0] = newPageIdLeft;
			nonLeafNode->pageNoArray[1] = newPageIdRight;

			//unpin
			bufMgr->unPinPage(file, newPageIdLeft, true);
			bufMgr->unPinPage(file, newPageIdRight, true);
			bufMgr->unPinPage(file, pageId, true);
			return;
		}

		//call recursive function
		T newChildValue;
		PageId newChildPageId = 0;
		insertEntryRecursive<T, T1, T2>
				(ridKeyPair,
				 nonLeafNode->pageNoArray[pos],
				 nonLeafNode->level == 1,
				 ARRAYMAX1,
				 ARRAYMAX2,
				 newChildValue,
				 newChildPageId);

		//check if child split
		if(newChildPageId != 0){
			//if child split

			//find last entry
			int last;
			for (last = 0; last < ARRAYMAX2; last++)
				if(nonLeafNode->pageNoArray[last] == 0)
					break;

			//check if full
			if(last < ARRAYMAX2){
				//not full
				for(int i=last; i>pos; i--){
					nonLeafNode->keyArray[i] = nonLeafNode->keyArray[i - 1];
					nonLeafNode->pageNoArray[i + 1] = nonLeafNode->pageNoArray[i];
				}
				nonLeafNode->keyArray[pos] = newChildValue;
				nonLeafNode->pageNoArray[pos + 1] = newChildPageId;
			}
			else{
				//full, need split
				Page* newPage;
				bufMgr->allocPage(file, newPageId, newPage);
				T2* newNonLeafNode = (T2*) newPage;

				//tmp array
				int tmpKeyArray[ARRAYMAX2 + 1];
				PageId tmpPageIdArray[ARRAYMAX2 + 2];

				//copy to tmp
				for(int i=0; i < ARRAYMAX2; i++) {
					tmpPageIdArray[i] = nonLeafNode->pageNoArray[i];
					tmpKeyArray[i] = nonLeafNode->keyArray[i];
					nonLeafNode->pageNoArray[i] = 0;
				}
				tmpPageIdArray[ARRAYMAX2 + 1] = nonLeafNode->pageNoArray[ARRAYMAX2 + 1];
				nonLeafNode->pageNoArray[ARRAYMAX2 + 1] = 0;

				for(int i= ARRAYMAX2; i > pos; i--){
					tmpKeyArray[i] = tmpKeyArray[i-1];
					tmpPageIdArray[i+1] = tmpPageIdArray[i];
				}
				tmpKeyArray[pos] = newChildValue;
				tmpPageIdArray[pos+1] = newChildPageId;

				//clear old and new page
				for(int i=0; i < ARRAYMAX2 + 1; i++){
					nonLeafNode->pageNoArray[i] = 0;
					newNonLeafNode->pageNoArray[i] = 0;
				}

				//copy back
				for(int i=0; i< (ARRAYMAX2 + 1) / 2; i++ ){
					nonLeafNode->keyArray[i] = tmpKeyArray[i];
					nonLeafNode->pageNoArray[i] = tmpPageIdArray[i];
				}
				nonLeafNode->pageNoArray[(ARRAYMAX2 + 1) / 2] = tmpPageIdArray[(ARRAYMAX2 + 1) / 2];

				for(int i= (ARRAYMAX2 + 1) / 2 + 1; i < ARRAYMAX2 + 1; i++){
					newNonLeafNode->keyArray[i] = tmpKeyArray[i];
					newNonLeafNode->pageNoArray[i] = tmpPageIdArray[i];
				}
				newNonLeafNode->pageNoArray[ARRAYMAX2 + 1] = tmpPageIdArray[ARRAYMAX2 + 1];

				//level
				newNonLeafNode->level = nonLeafNode->level;

				//push up
				newValue = tmpKeyArray[(ARRAYMAX2 + 1) / 2];

				//unpin
				bufMgr->unPinPage(file, newPageId, true);

			}
		}

		//unpin
		bufMgr->unPinPage(file, pageId, true);
	}
}

template<class T, class T1>
void BTreeIndex::handleNewRoot(T newValue, PageId newPageId, int ARRAYMAX){
	std::cout<<"handle new root"<<std::endl; //TODO
	PageId newRootPageId;
	Page *newRootPage;
	bufMgr->allocPage(file, newRootPageId, newRootPage);

	T1 *newRootNonLeafNode = (T1 *) newRootPage;
	for(int i=0;i<ARRAYMAX + 1; i++) newRootNonLeafNode->pageNoArray[i] = 0;
	newRootNonLeafNode->keyArray[0] = newValue;
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
		RIDKeyPair<char[STRINGSIZE] > ridKeyPair;

	}
	/*
		case STRING:{
			RIDKeyPair<char[STRINGSIZE] > ridKeyPair;
			strncpy(ridKeyPair.key, (char*) key, STRINGSIZE);
			ridKeyPair.rid = rid;
			char newValue[STRINGSIZE];
			PageId newPageId = 0;
//			insertEntryRecursive(ridKeyPair, rootPageNum, 0, newValue, newPageId);
			//TODO: root split


		}
		*/

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
		startScanHelper<int, NonLeafNodeInt>(*((int*) lowValParm), *((int *) highValParm));
	} else if (attributeType == DOUBLE) {
		this->lowValDouble = *((double*) lowValParm);
		this->highValDouble = *((double *) highValParm);
		startScanHelper<double , NonLeafNodeDouble>(*((double *) lowValParm), *((double *) highValParm));
	} else {
		//TODO
	}

}

template<class T, class T1>
void BTreeIndex::startScanHelper(T lowValParm,
								 T highValParm)
{
	if(lowValParm >  highValParm)
		throw BadScanrangeException();

	//find first one
	currentPageNum = rootPageNum;	std::cout<<"start scan read root page:" << currentPageNum<<std::endl;//TODO
	bufMgr->readPage(file, currentPageNum, currentPageData);
	T1* nonLeafNode = (T1*) currentPageData; std::cout<<"current page level:" << nonLeafNode->level<<std::endl;//TODO

	while(nonLeafNode->level != 1) {//std::cout<<"1";//TODO
		PageId nextPageId = nonLeafNode->pageNoArray[0];
		bufMgr->readPage(file, nextPageId, currentPageData);
		bufMgr->unPinPage(file, currentPageNum, false);
		currentPageNum = nextPageId;
		nonLeafNode = (T1*) currentPageData;
	}
	PageId nextPageId = nonLeafNode->pageNoArray[0];
	bufMgr->readPage(file, nextPageId, currentPageData);
	bufMgr->unPinPage(file, currentPageNum, false);
	currentPageNum = nextPageId;

	nextEntry = 0;

	std::cout<< "start scan complete" << std::endl;
}


// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

const void BTreeIndex::scanNext(RecordId& outRid) 
{
	if(scanExecuting == false){
		throw ScanNotInitializedException();
	}

	if(attributeType == INTEGER){
		scanNextHelper<int, LeafNodeInt> (outRid, lowValInt, highValInt, INTARRAYLEAFSIZE);
	} else if (attributeType == DOUBLE){
		scanNextHelper<double, LeafNodeDouble> (outRid, lowValDouble, highValDouble, DOUBLEARRAYLEAFSIZE);
	}

}

template <class T, class T1>
void BTreeIndex::scanNextHelper(RecordId &outRid, T lowVal, T highVal, int ARRAYMAX)
{
	T1* leafNode;
	while(1){
		leafNode = (T1*) currentPageData;
		if(leafNode->ridArray[nextEntry].page_number == 0
		   || nextEntry == ARRAYMAX) {

			PageId nextPageNum = leafNode->rightSibPageNo;
			if(nextPageNum == 0){
				std::cout<<"scan finish"<<std::endl; //TODO:delete
				throw IndexScanCompletedException();
			}

			bufMgr->unPinPage(file, currentPageNum, false);
			currentPageNum = nextPageNum;

			bufMgr->readPage(file, currentPageNum, currentPageData);

			nextEntry = 0;
			continue;
		}

		if((lowOp==GT && leafNode->keyArray[nextEntry] <= lowVal)
		   || (lowOp==GTE && leafNode->keyArray[nextEntry] < lowVal)
		   || (highOp==LT && leafNode->keyArray[nextEntry] >= highVal)
		   || (highOp==LTE && leafNode->keyArray[nextEntry] > highVal))
		{
			nextEntry++;
			continue;
		}

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
	if(scanExecuting == false){
		throw ScanNotInitializedException();
	}else {
		scanExecuting = false;
		bufMgr->unPinPage(file, currentPageNum, false);
		bufMgr->flushFile(file);
	}
}

}
