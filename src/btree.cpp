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
		indexMetaInfo->rootPageNo = rootPageNum;
		indexMetaInfo->attrByteOffset = attrByteOffset;
		indexMetaInfo->attrType = attrType;
		strcpy(indexMetaInfo->relationName, relationName.c_str());
		bufMgrIn->unPinPage(file, headerPageNum, true);

		//root page
		bufMgrIn->allocPage(file, rootPageNum, page);
		switch (attrType) {
			case INTEGER: {
				NonLeafNodeInt *nonLeafNodeInt = (NonLeafNodeInt *) page;
				nonLeafNodeInt->level = 1;
				break;
			}
			case DOUBLE: {
				NonLeafNodeDouble *nonLeafNodeDouble = (NonLeafNodeDouble *) page;
				nonLeafNodeDouble->level = 1;
				break;
			}
			case STRING: {
				NonLeafNodeString *nonLeafNodeString = (NonLeafNodeString *) page;
				nonLeafNodeString->level = 1;
				break;
			}
			default:
				break;
		}
		bufMgrIn->unPinPage(file, rootPageNum, true);

		//insert index
		FileScan fscan(relationName, bufMgr);
		RecordId scanRid;
		while(1)
		{
			fscan.scanNext(scanRid);
			insertEntry((void *) fscan.getRecord().c_str() + attrByteOffset, scanRid);	//TODO:?
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
		   || !strcmp(indexMetaInfo->relationName, relationName.c_str()))
			throw BadIndexInfoException("constructor parameters do not match exist index file");

		this->rootPageNum = indexMetaInfo->rootPageNo;
		bufMgrIn->unPinPage(file, headerPageNum, false);
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

void BTreeIndex::insertEntryRecursive(RIDKeyPair<int > ridKeyPair,
									  PageId pageId,
									  bool isLeaf,
									  int& newValue,
									  PageId& newPageId)
{
//	std::cout<<"try to insert key: "<< ridKeyPair.key << std::endl; //TODO:delete

	if(isLeaf){
		//read the page
		Page* page;
		bufMgr->readPage(file, pageId, page);
		LeafNodeInt* leafNodeInt = (LeafNodeInt*) page;

		//find position
		int pos = 0;
		while(ridKeyPair.key > leafNodeInt->keyArray[pos]
			  && leafNodeInt->ridArray[pos].page_number != 0
			  && pos < INTARRAYLEAFSIZE )
			pos++;

		//find last entry
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

//			std::cout<<"insert key: "<< ridKeyPair.key <<" at page :" << pageId << std::endl; //TODO:delete
		} else {
			//full
			Page* newPage;
			bufMgr->allocPage(file, newPageId, newPage);
			LeafNodeInt* newLeafNodeInt = (LeafNodeInt*) newPage;

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
				leafNodeInt->keyArray[i] = tmpKeyArray[i];
				leafNodeInt->ridArray[i] = tmpRidArray[i];
			}
			for(int i=(INTARRAYLEAFSIZE + 1)/2; i<INTARRAYLEAFSIZE + 1; i++){
				newLeafNodeInt->keyArray[i - (INTARRAYLEAFSIZE + 1)/2] = tmpKeyArray[i];
				newLeafNodeInt->ridArray[i - (INTARRAYLEAFSIZE + 1)/2] = tmpRidArray[i];
			}

			//link leaf node
			newLeafNodeInt->rightSibPageNo = leafNodeInt->rightSibPageNo;
			leafNodeInt->rightSibPageNo = newPageId;

			//push up
			newValue = newLeafNodeInt->keyArray[0];

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
		NonLeafNodeInt* nonLeafNodeInt = (NonLeafNodeInt*) page;

		//find pageArray position
		int pos = 0;
		while(ridKeyPair.key >= nonLeafNodeInt->keyArray[pos]
			  && nonLeafNodeInt->pageNoArray[pos+1] != 0
			  && pos < INTARRAYLEAFSIZE )
			pos++;

		//index file is empty
		if(nonLeafNodeInt->pageNoArray[pos] == 0){
			std::cout<<"empty index"  << std::endl; //TODO:delete
			PageId newPageIdLeft, newPageIdRight;
			Page* pageLeft, *pageRight;
			bufMgr->allocPage(file, newPageIdLeft, pageLeft);
			bufMgr->allocPage(file, newPageIdRight, pageRight);

			LeafNodeInt* leafNodeIntLeft = (LeafNodeInt*) pageLeft;
			LeafNodeInt* leafNodeIntRight = (LeafNodeInt*) pageRight;
			leafNodeIntRight->keyArray[0] = ridKeyPair.key;
			leafNodeIntRight->ridArray[0] = ridKeyPair.rid;
			leafNodeIntLeft->rightSibPageNo = newPageIdRight;
			leafNodeIntRight->rightSibPageNo = 0;
			nonLeafNodeInt->keyArray[0] = ridKeyPair.key;
			nonLeafNodeInt->pageNoArray[0] = newPageIdLeft;
			nonLeafNodeInt->pageNoArray[1] = newPageIdRight;

			//unpin
			bufMgr->unPinPage(file, newPageIdLeft, true);
			bufMgr->unPinPage(file, newPageIdRight, true);
			bufMgr->unPinPage(file, pageId, true);
			return;
		}

		//call recursive function
		int newChildValue = 0;
		PageId newChildPageId = 0;
		insertEntryRecursive(ridKeyPair,
							 nonLeafNodeInt->pageNoArray[pos],
							 nonLeafNodeInt->level==1,
							 newChildValue,
							 newChildPageId);

		//check if child split
		if(newChildPageId != 0){
			//if child split

			//find last entry
			int last;
			for (last = 0; last < INTARRAYLEAFSIZE; last++)
				if(nonLeafNodeInt->pageNoArray[last] == 0)
					break;

			//check if full
			if(last < INTARRAYLEAFSIZE){
				//not full
				for(int i=last; i>pos; i--){
					nonLeafNodeInt->keyArray[i] = nonLeafNodeInt->keyArray[i-1];
					nonLeafNodeInt->pageNoArray[i+1] = nonLeafNodeInt->pageNoArray[i];
				}
				nonLeafNodeInt->keyArray[pos] = newChildValue;
				nonLeafNodeInt->pageNoArray[pos+1] = newChildPageId;
			}
			else{
				//full, need split
				Page* newPage;
				bufMgr->allocPage(file, newPageId, newPage);
				NonLeafNodeInt* newNonLeafNodeInt = (NonLeafNodeInt*) newPage;

				//tmp array
				int tmpKeyArray[INTARRAYLEAFSIZE+1];
				PageId tmpPageIdArray[INTARRAYLEAFSIZE +2];

				//copy to tmp
				for(int i=0; i<INTARRAYLEAFSIZE; i++) {
					tmpPageIdArray[i] = nonLeafNodeInt->pageNoArray[i];
					tmpKeyArray[i] = nonLeafNodeInt->keyArray[i];
					nonLeafNodeInt->pageNoArray[i] = 0;
				}
				tmpPageIdArray[INTARRAYLEAFSIZE + 1] = nonLeafNodeInt->pageNoArray[INTARRAYLEAFSIZE +1];
				nonLeafNodeInt->pageNoArray[INTARRAYLEAFSIZE + 1] = 0;

				for(int i=INTARRAYLEAFSIZE; i>pos; i--){
					tmpKeyArray[i] = tmpKeyArray[i-1];
					tmpPageIdArray[i+1] = tmpPageIdArray[i];
				}
				tmpKeyArray[pos] = newChildValue;
				tmpPageIdArray[pos+1] = newChildPageId;

				//copy back
				for(int i=0; i<(INTARRAYLEAFSIZE+1)/2;i++ ){
					nonLeafNodeInt->keyArray[i] = tmpKeyArray[i];
					nonLeafNodeInt->pageNoArray[i] = tmpPageIdArray[i];
				}
				nonLeafNodeInt->pageNoArray[(INTARRAYLEAFSIZE+1)/2] = tmpPageIdArray[(INTARRAYLEAFSIZE+1)/2];

				for(int i=(INTARRAYLEAFSIZE + 1)/2 + 1; i<INTARRAYLEAFSIZE +1; i++){
					newNonLeafNodeInt->keyArray[i] = tmpKeyArray[i];
					newNonLeafNodeInt->pageNoArray[i] = tmpPageIdArray[i];
				}
				newNonLeafNodeInt->pageNoArray[INTARRAYLEAFSIZE + 1] = tmpPageIdArray[INTARRAYLEAFSIZE +1];

				//level
				newNonLeafNodeInt->level = nonLeafNodeInt->level;

				//push up
				newValue = tmpKeyArray[(INTARRAYLEAFSIZE+1) / 2];

				//unpin
				bufMgr->unPinPage(file, newPageId, true);

			}
		}

		//unpin
		bufMgr->unPinPage(file, pageId, true);
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

const void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{

	switch (attributeType){
		case INTEGER:{
			RIDKeyPair<int > ridKeyPair;
			ridKeyPair.set(rid, *((int *) key));
			int newValue;
			PageId newPageId = 0;
			insertEntryRecursive(ridKeyPair, rootPageNum, 0, newValue, newPageId);

			//if root got split
			if(newPageId != 0){
				PageId newRootPageId;
				Page* newRootPage;
				bufMgr->allocPage(file, newRootPageId, newRootPage);

				Page *page;
				bufMgr->readPage(file, rootPageNum, page);

				NonLeafNodeInt* newRootNonLeafNodeInt = (NonLeafNodeInt*) page;
				newRootNonLeafNodeInt->keyArray[0] = newValue;
				newRootNonLeafNodeInt->pageNoArray[0] = rootPageNum;
				newRootNonLeafNodeInt->pageNoArray[1] = newPageId;
				newRootNonLeafNodeInt->level = 0;
				rootPageNum = newRootPageId;
				bufMgr->unPinPage(file, newRootPageId, true);
				bufMgr->unPinPage(file, rootPageNum, true);
			}
		}
		case DOUBLE:{
			RIDKeyPair<double > ridKeyPair;
			ridKeyPair.set(rid, *((double *) key));
			double newValue;
			PageId newPageId = 0;
//			insertEntryRecursive(ridKeyPair, rootPageNum, 0, newValue, newPageId);
			//TODO: root split

		}
		case STRING:{
			RIDKeyPair<char[STRINGSIZE] > ridKeyPair;
			strncpy(ridKeyPair.key, (char*) key, STRINGSIZE);
			ridKeyPair.rid = rid;
			char newValue[STRINGSIZE];
			PageId newPageId = 0;
//			insertEntryRecursive(ridKeyPair, rootPageNum, 0, newValue, newPageId);
			//TODO: root split


		}
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

	switch (attributeType) {
		case INTEGER:{
			this->lowOp = lowOpParm;
			this->highOp = highOpParm;
			this->lowValInt = *((int*) lowValParm);
			this->highValInt = *((int *) highValParm);
		}
	}

	//find first one
	currentPageNum = rootPageNum;
	bufMgr->readPage(file, currentPageNum, currentPageData);
	std::cout<< "start scan, open rootpage"<< std::endl;	//TODO:delete
	NonLeafNodeInt* nonLeafNodeInt = (NonLeafNodeInt*) currentPageData;

	while(nonLeafNodeInt->level != 1) {
		PageId nextPageId = nonLeafNodeInt->pageNoArray[0];
		std::cout<<"start scan, go to page: "<< nextPageId <<std::endl; //TODO:endl
		bufMgr->readPage(file, nextPageId, currentPageData);
		bufMgr->unPinPage(file, currentPageNum, false);
		currentPageNum = nextPageId;
		nonLeafNodeInt = (NonLeafNodeInt*) currentPageData;
	}
	PageId nextPageId = nonLeafNodeInt->pageNoArray[0];
	bufMgr->readPage(file, nextPageId, currentPageData);
	bufMgr->unPinPage(file, currentPageNum, false);
	currentPageNum = nextPageId;
	nonLeafNodeInt = (NonLeafNodeInt*) currentPageData;

	nextEntry = 0;

	std::cout<< "start scan complete" << std::endl
	<<"current page:"<<currentPageNum << std::endl;	//TODO:delete
}

// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

const void BTreeIndex::scanNext(RecordId& outRid) 
{
	switch (attributeType){
		case INTEGER:{
			LeafNodeInt* leafNodeInt;
			while(1){
				leafNodeInt = (LeafNodeInt*) currentPageData;
				if(leafNodeInt->ridArray[nextEntry].page_number == 0
						|| nextEntry == INTARRAYLEAFSIZE) {
//					std::cout<<"scan next try to unpin page"<<currentPageNum << std::endl; //TODO:delete

					PageId nextPageNum = leafNodeInt->rightSibPageNo;
					if(nextPageNum == 0){
						std::cout<<"scan finish"<<std::endl; //TODO:delete
						throw IndexScanCompletedException();
					}

					bufMgr->unPinPage(file, currentPageNum, false);
					currentPageNum = nextPageNum;

					bufMgr->readPage(file, currentPageNum, currentPageData);
//					std::cout<<"scan next try to read page"<<currentPageNum << " complete"<<std::endl; //TODO:delete

					nextEntry = 0;
					continue;
				}

				if((lowOp==GT && leafNodeInt->keyArray[nextEntry] <= lowValInt)
						|| (lowOp==GTE && leafNodeInt->keyArray[nextEntry] < lowValInt)
						|| (highOp==LT && leafNodeInt->keyArray[nextEntry] >= highValInt)
						|| (highOp==LTE && leafNodeInt->keyArray[nextEntry] > highValInt))
				{
//					std::cout<< "not match: "<< leafNodeInt->keyArray[nextEntry] <<std::endl; //TODO:delete
					nextEntry++;
					continue;
				}

//				std::cout<< "match: "<< leafNodeInt->keyArray[nextEntry] <<std::endl; //TODO:delete
				outRid = leafNodeInt->ridArray[nextEntry];
				nextEntry++;
				return ;
			}
		}
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//
const void BTreeIndex::endScan() 
{
	scanExecuting = false;
	bufMgr->unPinPage(file, currentPageNum, false);
	bufMgr->flushFile(file);

}

}
