/****************************************************************************************

   Copyright (C) 2015 Autodesk, Inc.
   All rights reserved.

   Use of this software is subject to the terms of the Autodesk license agreement
   provided at the time of installation or download, or which otherwise accompanies
   this software in either electronic or hard copy form.

****************************************************************************************/
#ifndef _COMMON_H
#define _COMMON_H

#include <fbxsdk.h>

void initializeSdkObjects(FbxManager*& pManager, FbxScene*& pScene, std::string);
void destroySdkObjects(FbxManager* pManager, bool pExitStatus);

bool saveScene(FbxManager* pManager, FbxDocument* pScene, const char* pFilename, int pFileFormat=-1, bool pEmbedMedia=false);
bool loadScene(FbxManager* pManager, FbxDocument* pScene, const char* pFilename);

#endif // #ifndef _COMMON_H


