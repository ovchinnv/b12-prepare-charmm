% MATLAB build routines
%addpath('/home/taly/scripts/matlab/build');
addpath('util/matlab');
% check if this is octave :
qoct=exist('OCTAVE_VERSION','builtin');
%
name='1hzh';
qpsfgen=1; % for patch formatting

% segment ids for the protein chains : (optional)
chain_segid={ 'H', 'K', 'L' 'M' ; ...
              'HC1 ', 'HC2 ','LC1 ', 'LC2 ' };
% note that we need 4-character segids for correct reading of resulting PDB files in CHARMM

ch2seg= @(x) char(chain_segid(1+find(ismember(chain_segid(:), x)))); % get segment ID from chainID
ch2segt= @(x) strtrim(char(chain_segid(1+find(ismember(chain_segid(:), x))))); % get segment ID from chainID; trimmed to removed leading/trailing spaces

pdbfile=[name,'.pdb']; % in lieu of a "structure" file

if (~exist('qpdb')) ; qpdb=0 ; end

if (~qpdb) % read pdb file
% global molecule;
 disp(['==>Reading pdb structure from file ',pdbfile,' ...']);
 if (qoct); qoctpdb=1; elseif (~exist('qoctpdb')) ; qoctpdb=0 ; end
 if (qoctpdb)
  molecule=readpdb(pdbfile,1);% 1 for verbosity
 else
  molecule=pdbread(pdbfile);
 end
 qpdb=1; %
end
pdb=molecule.Model.Atom;

% manipulate PDB format and write out pdb files for charmm
anum= [pdb.AtomSerNo]';
aname={pdb.AtomName}' ;
rname={pdb.resName}' ;
segid={pdb.segID}' ;
resid=[pdb.resSeq]' ;
insertion=char({pdb.iCode});
occu =[pdb.occupancy]' ;
bet  =[pdb.tempFactor]' ;
xpdb =[pdb.X]';
ypdb =[pdb.Y]';
zpdb =[pdb.Z]';
chainid=[pdb.chainID];
element=[pdb.element];
natom=length(pdb);
% mark all chains (replacement of explicit code)
for ch = chain_segid(1,:)
 cch=char(ch)
 mcmd=['chain',cch,'=ismember(chainid,''',cch,'''); segid2=ch2seg(''',cch,''');for i=find(chain',cch,'); pdb(i).segID=segid2; end']
 eval(mcmd);
end
segid={pdb.segID}' ; % segid possibly redefined per above
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes to match charmm terminology
pdb=fixcharmm(pdb);
% create new model
mol2.Model.Atom=pdb;
% write out pdb files
if (1)
pdbout(mol2,[ch2segt('H'),'.pdb'],xpdb,ypdb,zpdb,[],[],find(chainH));
pdbout(mol2,[ch2segt('L'),'.pdb'],xpdb,ypdb,zpdb,[],[],find(chainL));
pdbout(mol2,[ch2segt('K'),'.pdb'],xpdb,ypdb,zpdb,[],[],find(chainK));
pdbout(mol2,[ch2segt('M'),'.pdb'],xpdb,ypdb,zpdb,[],[],find(chainM));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for possible disulfides via distance criterion
% optional selections (otherwise, all indices)
%seldisu1=ismember(strtrim(segid),'HC1')|ismember(strtrim(segid),'LC1')|ismember(strtrim(segid),'HC2')|ismember(strtrim(segid),'LC2');
%seldisu2=seldisu1;
disu ; % needs segid, resid, aname, rname to be defined above
%%%%%%% now process heterogeneous atom entries %%%%%%%%%%%%
pdbh=molecule.Model.HeterogenAtom;
anumh= [pdbh.AtomSerNo]';
anameh={pdbh.AtomName}' ;
rnameh={pdbh.resName}' ;
segidh={pdbh.segID}' ;
residh=[pdbh.resSeq]' ;
insertionh=char({pdbh.iCode});
occuh =[pdbh.occupancy]' ;
beth  =[pdbh.tempFactor]' ;
xpdbh =[pdbh.X]';
ypdbh =[pdbh.Y]';
zpdbh =[pdbh.Z]';
chainidh=[pdbh.chainID];
elementh=[pdbh.element];
%
% water
water=ismember(rnameh,'HOH');
resnum=1; % renumber residues
for i=find(water)'
 pdbh(i).resName='TP3';
 pdbh(i).AtomNameStruct.chemSymbol='O';
 pdbh(i).AtomNameStruct.remoteInd='H';
 pdbh(i).AtomNameStruct.branch='2';
 pdbh(i).segID='XWAT';
 pdbh(i).element='';
 pdbh(i).resSeq=resnum;%sprintf('%5d',resnum);
 resnum=resnum+1; % renumber residues
end

mol3.Model.Atom=pdbh;
% write out pdb files
pdbout(mol3,'xwat.pdb',xpdbh,ypdbh,zpdbh,[],[],find(water));

% rename NAG atoms to match CHARMM carbohydrate force field
%
nag=ismember(rnameh,'NAG');
nagc8=nag & ismember(anameh,'C8') ;
nago7=nag & ismember(anameh,'O7') ;
nagc7=nag & ismember(anameh,'C7') ;
nagn2=nag & ismember(anameh,'N2') ;

for i=find(nagc8)'
 pdbh(i).AtomNameStruct.chemSymbol='C';
 pdbh(i).AtomNameStruct.remoteInd='T';
end

for i=find(nago7)'
 pdbh(i).AtomNameStruct.chemSymbol='O';
 pdbh(i).AtomNameStruct.remoteInd='';
end

for i=find(nagc7)'
 pdbh(i).AtomNameStruct.chemSymbol='C';
 pdbh(i).AtomNameStruct.remoteInd='';
end

for i=find(nagn2)'
 pdbh(i).AtomNameStruct.chemSymbol='N';
 pdbh(i).AtomNameStruct.remoteInd='';
end

mol3.Model.Atom=pdbh;

% process glycans
glyco;
% create pdbs for missing loops
missing;
%
% write pdbs with missing coordinates included
%
pdbs={'HC2.pdb', 'HC2-MISSING-1.pdb', 'HC2-MISSING-2.pdb'};%, 'G120-MISSING-4.pdb'};
molout=combine_pdbs(pdbs);
pdbwrite('HC2-ALL.pdb', molout); system('echo END >> HC2-ALL.pdb');
%
