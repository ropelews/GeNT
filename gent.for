      Program GroupEntropy
C
      Implicit None
C
      Integer   NSP, MaxSqP, TotSqP, n20, n24, n25, StdIn, StdOut,
     +          DbgUnt, SqUnit, EnUnit, GpUnit, MxSymP, MxGrpP,
     +          NSPSqP
C
      Parameter  ( NSP = 600,   MaxSqP = 1000, TotSqP = NSP * MaxSqP,
     +             MxGrpP =25,  n24 = 24,    n20 = 20,   MxSymP = 28,
     +             StdIn = 5,   StdOut = 6,  SqUnit = 7, EnUnit = 8,
     +             GpUnit = 9,  DbgUnt = 10, n25 = 25,
     +             NSPSqP = ( ( NSP * (NSP+1) ) / 2 )     )
C
      Integer         KeyBrd, Termnl, SeqFl, GrupFl, OutFl, DeBugF,
     +                AlnLen, DoGrup, NotID, ModAA, PFLen, WKnt, BKnt,
     +                LNoGap, NSqUsed, Unknown, n, LsTot, NS, NGrp1,
     +                ASize, Gap, MaxLet, GrpKnt, l, MatSiz, KMotif,
     +                KMPos
      Integer         LSN( NSP ), LsCum( NSP ), Seq( TotSqP ),
     +                GrupId( NSP ), InvGrp( NSP, MxGrpP ), AlfMap(20),
     +                Group1( NSP ), LGpNam( MxGrpP ), MapAA( 20 ),
     +                NGroup( MxGrpP ), Comp(0:MxSymP), MapNuc(4),
     +                NoGaps(MaxSqP), GapMap(MaxSqP), GComp(0:MxSymP),
     +                SeqUsed( NSP ), NMotif( MaxSqP ), IMotif(MaxSqP),
     +                LMotif( MaxSqP ), BMotif( MaxSqP )
      Real*4          MaxMat, MinMat, MPC, AlnEnt, AlnScr, MatAve,
     +                WAve, BAve, Wstd, Bstd, FZAv, FZStd, GZAv, GZStd,
     +                ZLimit
      Real*4          SeqWt( NSP ), Profil( 25, MaxSqP, 2, MxGrpP ),
     +                ProfileT(23, MaxSqP), Cross( MaxSqP, 2, MxGrpP ),
     +                Entropy( MaxSqP ), Qij( n20, n20 ), QSum( n20 ),
     +                SwProt(n20), CVPScr(TotSqP,2), CVPEnt(TotSqP,2),
     +                CVScor(NSP,2), CVEntr(NSP,2), TotCE( MxGrpP, 2 ),
     +                RawKnt(25,MaxSqP), BackF(25,MaxSqP), FScore(NSP),
     +                MGpScr(MxGrpP), thresh(MxGrpP,2), Within(NSPSQP),
     +                GrpAln( NSP, 0:MxGrpP ), Pair( NSP, NSP ),
     +                AlOthr( MaxSqP, 2, MxGrpP ), Between( NSPSQP ),
     +                BackFr(MxSymP), FZEnt(MaxSqP),  MaxGZ(MaxSqP),
     +                GZEnt(MaxSqP,MxGrpP), FrGap(MaxSqP),
     +                MaxVol(MaxSqP), MinVol(MaxSqP), AveVol(MaxSqP),
     +                SigVol( MaxSqP )
      Character*1     BestSq( MaxSqP, 2, MxGrpP ), BestSqT( MaxSqP ),
     *                sqprnt(25), MemeAA(20), MemNuc(4), Alfbet(20),
     +                UCPrnt(25), GCons(MaxSqP,2,MxGrpP),
     +                Consen(MaxSqP)
      Character*4     Bits
      Character*15    SqName( NSP )
      Character*24    PFile
      Character*30    GrpNam( MxGrpP )
      Character*60    MatDef
      Character*70    Title
      Character*80    Line
      Logical         Protin, DNA, empty, Split, DoPair, AsnGrp, Pass1,
     +                HaveGp, CVout, PlotOut, GDout, PSSMout
C
C
C AJR2020 temp vars for debug
      integer it, i, j




      Data MemeAA  / 'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l',
     +               'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y' /
      Data MapAA   / 1, 5, 4, 7, 14, 8, 9, 10, 12, 11, 13, 3, 15, 6,
     +               2, 16, 17, 20, 18, 19 /
      Data MemNuc  / 'a', 'c', 'g', 't'  /
      Data MapNuc  / 1, 2, 3, 4 /
C
      Data UCPrnt / 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
     +              'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
     +              'B', 'Z', 'X', '-', '.' /
C
C
C ***   Order of the amino acids for the array SwProt, the fractional
C ***   composition is derived from PIR release 56.
C
C    Array         One letter Amino Acid codes
C    Index
C    1 - 14      A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,  F,
C   15 - 20      P,  S,  T,  W,  Y,  V
C
      Data  SwProt / 0.074188, 0.050964, 0.045248, 0.052391, 0.016330,
     +               0.040148, 0.064263, 0.068571, 0.021808, 0.058507,
     +               0.094275, 0.060616, 0.023182, 0.040902, 0.048982,
     +               0.072478, 0.056846, 0.012443, 0.032713, 0.065144 /
C
C
    1 Format(//' ***   Classification and Evaluation of Sequences in',
     *                 ' Biochemical Groups  ***',
     1       //'       This program computes the cross entropy for',
     2                 ' groups of sequences',
     3        /'       that have been assigned to groups on the basis',
     4                 ' of biochemical,',
     5        /'       physiological, or other biological property.',
     6                 '  The sequence',
     7        /'       assignments are cross-validated, again by the',
     8                 ' cross entropy',
     9        /'       measure, to check for problems with the',
     A                 ' alignment or group',
     B        /'       assignment.  Sequences that were initially',
     C                 ' identified as',
     D        /'       "unclassified" are compared to all of the',
     E                 ' groups using position',
     F        /'       specific log-odds scores as described by',
     G                 ' Henikoff and Henikoff.',
     H        /'       Positions in the aligned sequences that are',
     I                 ' important for',
     J        /'       determining group membership are identified by',
     K                 ' having a high',
     L        /'       entropy for the entire alignment and a high',
     M                 ' entropy for one or',
     N        /'       more specific groups.', / )
    2 Format( A80 )
    4 Format(//' Computing Entropy for the entire alignment.')
    5 Format(//' Computing Entropy for Group: ', A30 )
    6 Format(//' Suming Other Group Cross-Validations.' )
    7 Format(//' Computing Entropy Distances between all pairs of',
     *         ' sequences.' )
    8 Format(//' Writing Cross-Validation Report files.' )
    9 Format(//' Writing Plot files of data to identify sequence',
     *         ' positons diagnostic of',
     1        /' membership in a specific group.')
   10 Format(//' Writing a file of data to use in assigning sequences',
     *         ' to groups.' )
   11 Format(//' Writing Position Specific Scoring Matrices for each',
     *         ' group to a file.' )
   12 Format(//' Writing distances between pairs of sequences to a',
     *         ' file.' )
   14 Format(//' Should the alignment scores used for interactive',
     *         ' classification of',
     1        /' sequences use alignment scores for all analyzed',
     2         ' sequence positions,',
     3        /' or use a limited set of positions that have a high',
     4         ' information',
     5        /' content for at least one of the defined groups. ',
     6         ' Answer A for ALL',
     7        /' sequence positions or L for a limited set of',
     8         ' sequence positions.', / )
   15 Format(//' What is the allowed lowest Z score for the',
     *         ' positions used to classify',
     1        /' currently unclassified sequences?', / )
   16 Format(//' Do you want to write a file of the new group',
     *         ' assignments?  Answer',
     1        /' (Y)es or (N)o.  The default is Yes.' )
   17 Format(//' Do you want to do interactive assignment of the so',
     *         ' far unclassified',
     1        /' sequences to the already defined groups or would you',
     2         ' prefer to do the',
     3        /' group assignment later and off-line.  Answer',
     4         ' (I)nteractive or (O)ff-line.',
     5        /' The default is Off-line.' )
   18 Format(//' Do you want the program to compute an entropy',
     *         ' distance between every pair of',
     1        /' sequences in the alignment?  Answer (Y)es or (N)o.',
     2         '  The default is No.' )
   19 Format(//' Select the results files you want written by entering',
     *         ' the numbers',
     1        /' corresponding to the kinds of results files you want',
     2         ' written.',
     3       //' 1.  Cross-validation files -- check groups for',
     4               ' nonconforming members.',
     5        /' 2.  Plot files -- data for family and group entropies',
     6               ' by sequence position.',
     7        /' 3.  GeneDoc files  --  highlight an alignment at high',
     8               ' entropy positions.',
     9        /' 4.  Profile file -- use to search a database or',
     A               ' identify group members.' )
   20 Format(//' *****   ERROR   *****',/
     *        /' Your data contains ',I4,' sequences of length ',I6,
     1         ' residues.',
     2        /' The GEnt program is currently configured to accept,',
     3         ' at most', I5,
     4        /' sequences of no more than ', I6,' residues.'/
     5        /' This can be corrected either by trimming the',
     6         ' N-terminal and C-terminal',
     7        /' ends of the sequences to remove positions containing',
     8         ' alignment gaps.',
     9        /' These positions would be ignored by the GEnt program',
     A         ' during its analysis',
     B        /' of your alignment.  Alternatively, send e-mail to:',
     C        /'           nicholas@psc.edu',
     D        /' requesting that the program be reconfigured accept',
     E         ' your alignment.',// )
   34 Format(//' Please enter the pseudocount multiplier, the',
     *         ' default value is 5.0',/)
C
C
      KeyBrd = StdIn
      Termnl = StdOut
      SeqFl = SqUnit
      GrupFl = GpUnit
      OutFl = EnUnit
      DeBugF = DbgUnt
      DoPair = .False.
      Pass1 = .True.
      HaveGp = .True.
      Write( Termnl, 1 )
C
      Open( Unit = OutFl, File = 'GroupClassifier.Out', Status = 'new',
     +      Access = 'Sequential', Form = 'formatted')
C     +      CarriageControl = 'list' )
C
C      Open( Unit = DeBugF, File = 'GroupClassifier.DeBug',
C     *      Access = 'Sequential', Form = 'formatted',
C     1      CarriageControl = 'list', Status = 'new' )
C
C**********************************************************************
C AJR IMPORTANT NOTES ABOUT THIS CODE AND VARIABLES. 
C NOT SURE THIS IS 100% ACCURATE AS IT WAS DOCUMENTED
C THROUGH TRIAL AND ERROR
C     SEQ = contains the sequences, linearized as a one-dimensional 
C           array, converted into an interger representation
C     LSTOT = Total length of all the sequences (e.g. the amount of 
C           actual sequence information stored in SEQ)
C     LSN = Array containing the length of each sequence N. Can be 
C           used to index into SEQ array. 
C     NS =  the number of sequences read in. (e.g. the amount of
C           useful data in array LSN)      
C     UCPRINT = the character table that converts the numeric sequence 
C           representation into the actual sequence.
C           to print out the sequences, use code like:
C                it=1
C                do 22 i=1,ns
C                   write(6,*) sqname(i),lsn(i),it
C                   write(6,*) ( ucprnt(seq(j)), j=it, it-1+lsn(i), 1)
C                   it=it+lsn(i)
C           22   continue

C           
C**********************************************************************

      Call GetSeq( Termnl, KeyBrd, SeqFl, OutFl, NSP, TotSqP, LsTot,
     +             NS, SEQ, LSN, LsCum, SqName, SeqWt, MxSymP, Comp,
     +             Gap, MaxLet, SqPrnt, n25, ModAA, Protin, DNA,
     +             DeBugF, AlnLen, ASize, BackFr, Unknown )
      If( NS .gt. NSP  .or.  AlnLen .gt. MaxSqP )    Then
         Write( Termnl, 20 )  NS, AlnLen, NSP, MaxSqP
         Close( Unit = OutFl, status = 'keep' )
         Stop ' Stopped for data too large.'
      Else
      EndIf
CAJR -- Fixed to raise this error if all sequences are not the same length'
      do 33 i=1,NS,1
        if (lsn(i) .ne. lsn(1)) then
           write(6,*) 
     +     'Sequence lengths are not equal -- are they aligned?'
           stop ' Stopped for sequences not being aligned.'
        endif
 33   continue
C
CAJR2020 This will print out the sequences read in:
C 
C      it=1
C      do 22 i=1,ns
C         write(6,*) sqname(i),lsn(i),it
C         write(6,*) ( seq(j), j=it, it-1+lsn(i), 1)
C         write(6,*) ( ucprnt(seq(j)), j=it, it-1+lsn(i), 1)
C         it=it+lsn(i)
C 22   continue

C      write(6,*) 'lsn=',lsn
C      write(6,*) 'ns=',ns
C      write(6,*) 'NSP=',nsp
C      write(6,*) 'seq=',Seq
C      write(6,*) 'sqprnt=',SqPrnt
      Call GetGrp( Termnl, KeyBrd, GrupFl, NS, NSP, SqName, MxGrpP,
     *             GrupID, InvGrp, GrpNam, LGpNam, NGroup, Title,
     1             PFile, PFLen, NotID, GrpKnt, OutFl, DeBugF, NSqUsed,
     2             SeqUsed )
C
C      Call RMotif( Termnl, Keybrd, MaxSqP, NMotif, IMotif, LMotif,      Public
C     *             BMotif, KMotif, KMPos )                              Public
C
      Call RmGaps( MaxSqP, TotSqP, NSP, NS, MxSymP, AlnLen, Gap,
     +             LNoGap, LsCum, Seq, GrupId, NoGaps, GapMap, NSqUsed,
     +             FrGap, BackFr, ASize, UnKnown, GComp, OutFl, MaxVol,
     +             MinVol, AveVol, SigVol, Consen, SqPrnt, Termnl,
     +             KeyBrd, GrupFl, SqName )
C
      Call GetQij( Termnl, KeyBrd, OutFl, Bits, Qij, QSum, n20, n24,
     *            MatAve, MaxMat, MinMat, MatSiz, MatDef, SqPrnt )
C
      Write( Termnl, 34 )
      Read( KeyBrd, * )   MPC
      If( MPC .le. 0.00001 )    MPC = 5.000
C
  200 Write( Termnl, 17 )
      Read( KeyBrd, 2 )  Line
      Call PakNam( Line, 80, empty )
      If( Line(1:1) .eq. 'I'  .or.  Line(1:1) .eq. 'i' )    Then
         AsnGrp = .True.
         Write( Termnl, 14 )
         Read( KeyBrd, 2 )  Line
         Call PakNam( Line, 80, empty )
         If( Line(1:1) .eq. 'L'  .or.  Line(1:1) .eq. 'l' )    Then
            Write( Termnl, 15 )
            Read( KeyBrd, * )  ZLimit
         Else
            ZLimit = -1.0
         EndIf
      Else
         AsnGrp = .False.
      EndIf
C
      Write( Termnl, 18 )
      Read( KeyBrd, 2 )  Line
      Call PakNam( Line, 80, empty )
      If( Line(1:1) .eq. 'Y'  .or.  Line(1:1) .eq. 'y' )    Then
         DoPair = .True.
      Else
         DoPair = .False.
      EndIf
C
      CVout = .False.
      PlotOut = .False.
      GDout = .False.
      PSSMout = .False.
      Write( Termnl, 19 )
      Read( KeyBrd, 2 )  Line
      If( Index( Line, '1' ) .gt. 0 )    CVout = .True.
      If( Index( Line, '2' ) .gt. 0 )    PlotOut= .True.
      If( Index( Line, '3' ) .gt. 0 )    GDout = .True.
      If( Index( Line, '4' ) .gt. 0 )    PSSMout = .True.
C
      Write( Termnl, 4 )
      Call EWhole( NSP, TotSqP, MaxSqP, n20, NS, Seq, LsCum, GrupID,
     *             ProfileT, Qij, QSum, Entropy, ASize, AlnLen, NoGaps,
     1             LNoGap, Protin, MPC, SwProt, BestSqT, MaxLet, SqPrnt,
     2             AlnEnt, AlnScr, FScore, DeBugF, NotID, Consen )
C
      Do 450 DoGrup = 1, GrpKnt, 1
         Do 410 l = 1, MaxSqP, 1
            AlOthr( l, 1, DoGrup ) = 0.0
            AlOthr( l, 2, DoGrup ) = 0.0
  410       continue
         If( DoGrup .eq. NotID )    GoTo 450
C
         Write( Termnl, 5 )  GrpNam( DoGrup )
         Call GrpEnt( NSP, TotSqP, MaxSqP, n20, NS, MxSymP, Seq,
     *                LsCum, GrupID, Qij, QSum, Profil, Cross, RawKnt,
     1                BackF, AlnLen, NoGaps, LNoGap, ASize, Protin,
     2                MPC, BestSq, GCons, Gap, MaxLet, SqPrnt, MxGrpP,
     3                NotID, DoGrup, MGpScr, GrpAln, DeBugF, TotCE )
C
         Do 430 n = 1, NGroup( DoGrup ), 1
            Group1( n ) = InvGrp( n, DoGrup )
  430       continue
         NGrp1 = NGroup( DoGrup )
C
         Call CRVal( NSP, TotSqP, MaxSqP, n20, MxSymP, SEQ, LsCum,
     *               Group1, Qij, QSum, CVPEnt, CVPScr, CVEntr,
     1               CVScor, RawKnt, BackF, AlnLen, NoGaps, LNoGap,
     2               ASize, MaxLet, Protin, MPC, NGrp1, SqPrnt, MxGrpP,
     3               thresh, DoGrup, unknown, BackFr, DeBugF )
C
  450    continue
C
C ***   Write files containing the group consensus (most frequent residue)
C ***   and the group identifier (most informative or highest entropy
C ***   residue) sequences.
C
C      Call WrCons( SeqFl, AlnLen, GrpKnt, PFLen, PFile, MaxSqP, MxGrpP,
C     *             Consen, GCons, BestSqT, BestSq, GrpNam )
C
C ***   Sum the cross entropies for all groups other than the group
C ***   currently being analyzed
C
      Write( Termnl, 6 )
      Call CVOthr( GrpKnt, LNoGap, NoGaps, Cross,AlOthr,MxGrpP,MaxSqP)
C
C ***   Compute the Z score at each position for the family entropy and
C ***   all of the group entropies.  Find the maximum Z score over all
C ***   groups at each position.
C
      Call ZScore( OutFl, MaxSqP, MxGrpP, NSP, TotSqP, Seq, LsCum,
     *             AlnLen, LGpNam, Entropy, Cross, GrpNam, SqPrnt,
     1             LNoGap, GrpKnt, InvGrp, SqName, GapMap, GZav,
     2             GZstd, FZav, FZstd, MaxGZ, FZEnt, GZEnt, NotId )
C
C ***   Compute the entropy (Information) distance between every pair
C ***   of sequences.  For every pair of sequences where both sequences
C ***   are members of a defined group, put the distances into either
C ***   the within group distribution or the between group distribution.
C
      If( DoPair )    Then
         Write( Termnl, 7 )
         Call TwoEnt( NSP, TotSqP, MaxSqP, n20, NS, MxSymP, Seq, LsCum,
     *                Qij, QSum, Pair, WithIn, Between, AlnLen, NoGaps,
     1                LNoGap, ASize, WKnt, BKnt, Protin, MPC, GrupID,
     2                NotID, Gap, MaxLet, BackFr, UnKnown, DeBugF,
     3                NSPSqP )
C
         Call AveVar( WithIn, WKnt, WAve, Wstd )
         Call AveVar( Between, BKnt, BAve, Bstd )
      Else
      EndIf
C
C ***   If the user has requested it do interactive assignment of
C ***   unclassified sequences to defined groups and cycle back
C ***   through the calculations rather than writing final files.
C
      If( AsnGrp )    Then
         Pass1 = .False.
         Call IAsnGp( Termnl, Title, GrpNam, GrpKnt, NGroup, InvGrp,
     *                SqName, Fscore, GrpAln, GrupId, NS, NSP,
     1                MxGrpP, MaxSqP, NSqUsed, SeqUsed, KeyBrd,
     2                Profil, ProfileT, GZEnt, MaxGZ, ZLimit,
     3                TotSqP, LNoGap, NoGaps, LsCum, Seq )
         Write( Termnl, 16 )
         Read( KeyBrd, 2 )  Line
         Call PakNam( Line, 80, empty )
         If( Line(1:1) .eq. 'N'  .or.  Line(1:1) .eq. 'n' )    Then
            HaveGp = .False.
         Else
            Call WrGpFl( Termnl, KeyBrd, GrupFl, NSP, MxGrpP, SqName,
     *                   GrpKnt, GrpNam, Title, PFile, InvGrp, NGroup)
            HaveGp = .True.
         EndIf
         GoTo 200
      Else
      EndIf
C
C ***   Write the Cross-validation reports, one file for each group,
C ***   named for the group with the extension .cvg
C
      If( CVout )    Then
         Write( Termnl, 8 )
         Call WritCV( GrupFl, NotID, GrpKnt, GrpNam, LGpNam, Title,
     *                NGroup, InvGrp, NoGaps, LNoGap, AlnEnt, TotCE,
     1                CVEntr, AlnScr, MGpScr, CVScor, Entropy, Cross,
     2                LsCum, CVPEnt, BestSq, BestSqT, NSP, MxGrpP,
     3                MaxSqP, TotSqP )
      EndIf
C
C ***   Write the Plot report files - data to visually examine which
C ***   sequence positions best discriminate a group from the others.
C ***   Assumes that the highest numbered group is the sequences that
C ***   are not assigned to a defined group.
C
      If( PlotOut )    Then
         Write( Termnl, 9 )
         Call WrPlot( GrupFl, GrpKnt, NotID, LGpNam, GrpNam, Title,
     *                LNoGap, NoGaps, BestSqT, Entropy, BestSq,
     1                Cross, AlOthr, NSP, MxGrpP, MaxSqP )
      EndIf
C
C ***   Write the data to assign ungrouped sequences to one of the user
C ***   predefined groups.
C
      Write( Termnl, 10 )
      Call WAsnGp( GrupFl, Title, GrpNam, GrpKnt, NGroup, InvGrp,
     *             SqName, Fscore, GrpAln, CVScor, GrupId, NS,
     1             NSP, MxGrpP )
C
C ***   Write the profiles or Possition Specific Scoring Matrices to
C ***   a file
C
      If( PSSMout )    Then
         Write( Termnl, 11 )
         Call WrPSSM( Termnl, GrupFl, ProfileT, Profil, MaxSqP,
     *                MxGrpP, AlfMap, MapAA, MapNuc, NoGaps, MemeAA,
     1                MemNuc, Alfbet, GrpNam, Title, Protin, PFile,
     2                ASize, LNoGap, thresh, GrpKnt, PFLen, AlnLen,
     3                MxSymP, Comp )
      EndIf
C
C ***   Write the Within and Between distance distributions to files
C ***   and create a Phylip formatted file of symmetrical cross-
C ***   entropy distances between each pair of sequences.
C
      If( DoPair )    Then
         Write( Termnl, 12 )
         Call TwoOut( Termnl, GrupFl, NSPSqP, NSP, WKnt, BKnt, PFLen,
     *                PFile, NSqUsed, SeqUsed, SqName, WAve, BAve,
     1                Wstd, Bstd, Within, Between, Pair )
      Else
      EndIf
C
C ***   Write files to import into GeneDoc for highlighting a multiple
C ***   sequence alignment
C
      If( GDout )    Then
         Call GeneDoc( GrupFl, MaxSqP, MxGrpP, NSP, TotSqP,
     *                 Gap, Seq, LsCum, AlnLen, LGpNam, Entropy,
     1                 Cross, GrpNam, SqPrnt, LNoGap, GrpKnt,
     2                 InvGrp, SqName, GapMap, NotID, GZEnt, FZEnt )
      EndIf
C
C ***  If the last interactive assignment of unclassified sequences to
C ***  groups was not saved, write it to a file now.
C
      If( .not. HaveGp )    Then
         Call WrGpFl( Termnl, KeyBrd, GrupFl, NSP, MxGrpP, SqName,
     *                GrpKnt, GrpNam, Title, PFile, InvGrp, NGroup )
         HaveGp = .True.
      Else
      EndIf
C
      Close( Unit = OutFl, Status = 'keep' )
C
C ***   Write volume and conservation analysis for each column of the
C ***   alignment.
C
C      Call WrVol( MaxSqP, TotSqP, NSP, NS, MxSymP, AlnLen, Gap, LNoGap,
C     *            LsCum, Seq, GrupId, NoGaps, GapMap, NSqUsed, FrGap,
C     1            BackFr, ASize, UnKnown, NMotif, Entropy, Consen,
C     3            BestSqT, sqprnt, MaxVol, MinVol, AveVol, SigVol )
C
      Stop 'Group Cross Entropies Calculated.'
      End
C
      Subroutine CVOthr( GrpKnt, LNoGap, NoGaps, Cross, AlOthr, MxGrpP,
     *                   MaxSqP )
C
C ***   Cross-Validate Others group
C
      Implicit None
C
      Integer     GrpKnt, LNoGap, MxGrpP, MaxSqP
      Integer     NoGaps(MaxSqP)
      Real*4      Cross( MaxSqP, 2, MxGrpP ), AlOthr( MaxSqP,2,MxGrpP)
C
      Integer    DoGrup, n, la, l
C
      Do 540 DoGrup = 1, GrpKnt, 1
         Do 530 n = 1, GrpKnt, 1
            If( n .eq. DoGrup )    GoTo 530
            Do 520 la = 1, LNoGap, 1
               l = NoGaps(la)
               AlOthr(l,1,DoGrup) = AlOthr(l,1,DoGrup) + Cross(l,1,n)
               AlOthr(l,2,DoGrup) = AlOthr(l,2,DoGrup) + Cross(l,2,n)
  520          continue
  530       continue
  540    continue
C
      Return
      End
C
      Subroutine WritCV( GrupFl, NotID, GrpKnt, GrpNam, LGpNam, Title,
     *                   NGroup, InvGrp, NoGaps, LNoGap, AlnEnt, TotCE,
     1                   CVEntr, AlnScr, MGpScr, CVScor, Entropy, Cross,
     2                   LsCum, CVPEnt, BestSq, BestSqT, NSP, MxGrpP,
     3                   MaxSqP, TotSqP )
C
C  ***   Write the Cross-Validation results to file for each group
C
      Implicit None
C
C ***   Passed Variable Declarations
C
      Integer        GrupFl, NSP, MxGrpP, MaxSqP, NotID, LNoGap, GrpKnt,
     *               TotSqP
      Integer        LGpNam( MxGrpP ), NoGaps(MaxSqP), NGroup(MxGrpP),
     *               InvGrp( NSP, MxGrpP ), LsCum( NSP )
      Real*4         AlnEnt, AlnScr
      Real*4         Entropy( MaxSqP ), Cross( MaxSqP, 2, MxGrpP ),
     *               MGpScr(MxGrpP), CVScor(NSP,2), CVEntr(NSP,2),
     1               TotCE( MxGrpP, 2 ), CVPEnt( TotSqP, 2 )
      Character*1    BestSq( MaxSqP, 2, MxGrpP ), BestSqT( MaxSqP )
      Character*30   GrpNam( MxGrpP )
      Character*70   Title
C
C ***   Local Variable Declarations
C
C                    LGrupP must be as large as global parameter MXGrpP
      Integer        LGrupP
      Parameter   (  LGrupP = 25 )
      Integer        l, la, DoGrup, NBlock, nb, i, ib, ie, RecLen
      Real*4         TmpEnt( LGrupP )         
      Character*1    TmpCh( LGrupP )
      Character*4    ch4
      Character*34   CVFile
      Character*44   VFmt
C
      Data  VFmt
     *    / '('' '',I3,'' '',I5,'' '',A1,xxxx('' '',I4,'' '',A1:)) ' /
C            1 23 45678 9_ 12345 67 89_1234567 89 _1234 56 789_1234
C
    2 Format(/' ', A70, / )
    3 Format(' Group: ', I4, ' ... ', A30 )
   11 Format( ' ', I6, F10.3, ' ', A1, '  ', A1, 3F10.3, 8( F10.3:) )
   12 Format(/' Field 1: Alignment Position or Index',
     *    /' Field 2: Entropy at Position for all sequences',
     1    /' Field 3: Consensus Sequence for the entire alignment',
     2    /' Field 4: Consensus Sequence for the group alignment',
     3    /' Field 5: Cross Entropy at Position for group (total)',
     4    /' Field 6: Cross Entropy at Position for group ( p(i))',
     5    /' Field 7: Cross Entropy at Position for not-group (q(i))',
     6    /' Fields 8 -> 17: Total cross entropy for the group at',
     7                       ' the position a',
     8    /'                 sequence removed from the group', / )
   13 Format(/' Field 1: Entropy for all aligned sequences',
     *    /' Field 2: Cross Entropy for group (total) alignment',
     1    /' Field 3: Cross Entropy component for group ( p(i))',
     2    /' Field 4: Cross Entropy component for not-group (q(i))',
     3    /' Fields 5 -> 14: Total cross entropy for the group',
     4                       ' with a sequence',
     4    /'                 removed from the group', / )
   14 Format( ' ', F10.3, ' .. ', 3F10.3, ' ... ', 8( F10.3:) )
   15 Format(/' Field 1: Best possible score for all sequences',
     *    /' Field 2: Best possible score for group alignment',
     1    /' Fields 3 -> 12: Cross validation score for the sequence',
     2                       ' against the group', / )
   16 Format( ' ', F10.3, ' .. ', F10.3, ' ... ', 8( F10.3:) )
   17 Format( I4 )
   18 Format('//      Entropy for the alignment and Cross Entropies',
     *       ' for all groups',
     1      /'//       followed by the most informative residue for',
     2       ' group membership.', /'// ', A70 )
   19 Format('// Data Format:  <',A44,'>' )
C
      Do 900 DoGrup = 1, GrpKnt, 1
         If( DoGrup .eq. NotID )    GoTo 900
         CVFile(1:30) = GrpNam( DoGrup )
         CVFile(31:34) = '    '
         CVFile(LGpNam(DoGrup)+1:LGpNam(DoGrup)+4) = '.cvg'
         Open( Unit = GrupFl, File = CVFile, Status = 'new',
     *         Access = 'Sequential', form = 'Formatted')
c     *         CarriageControl = 'List' )
         Write( GrupFl, 2 )  Title
         Write( GrupFl, 3 )  DoGrup, GrpNam( DoGrup )
         NBlock = ( NGroup( DoGrup ) + 7 ) / 8
         ib = -7
         Write( GrupFl, 13 )
         Do 700 nb = 1, NBlock, 1
            ib = ib + 8
            ie = ib + 7
            If( ie .gt. NGroup(DoGrup) )    ie = NGroup(DoGrup)
            Write( GrupFl, 14 ) AlnEnt,
     *                          TotCE(DoGrup,1) + TotCE(DoGrup,2),
     1                          TotCE(DoGrup,1), TotCE(DoGrup,2),
     2                         (( CVEntr(i,1) + CVEntr(i,2)),
     3                            i = ib, ie, 1 )
  700       continue
         ib = -7
         Write( GrupFl, 15 )
         Do 750 nb = 1, NBlock, 1
            ib = ib + 8
            ie = ib + 7
            If( ie .gt. NGroup(DoGrup) )    ie = NGroup(DoGrup)
            Write( GrupFl, 16 ) AlnScr, MGpScr(DoGrup),
     1                         ( CVScor(i,1), i = ib, ie, 1 )
  750       continue
         ib = -7
         Write( GrupFl, 12 )
         Do 850 nb = 1, NBlock, 1
            ib = ib + 8
            ie = ib + 7
            If( ie .gt. NGroup(DoGrup) )    ie = NGroup(DoGrup)
            Do 800 la = 1, LNoGap, 1
               l = NoGaps(la)
               Write( GrupFl, 11 ) l, Entropy(l), BestSqT(l),
     *                             BestSq(l,1,DoGrup),
     1                             Cross(l,1,DoGrup)+Cross(l,2,DoGrup),
     2                             Cross(l,1,DoGrup),Cross(l,2,DoGrup),
     3                         (( CVPEnt(LsCum(InvGrp(i,DoGrup))+l,1)
     4                          + CVPEnt(LsCum(InvGrp(i,DoGrup))+l,2)),
     5                            i = ib, ie, 1 )
  800          continue
  850       continue
         Close( Unit = GrupFl, Status = 'keep' )
  900    continue
C
C      RecLen = 13 + (( GrpKnt-1 ) * 7 )
C      Open( Unit = GrupFl, File = 'AllGroups.cvg', Status = 'new',
C     *      Access = 'Sequential', form = 'Formatted',
C     *      CarriageControl = 'List', RECL = RecLen )
C      Write( ch4, 17 )  GrpKnt-1
C      VFmt(23:26) = ch4
C      Write( GrupFl, 19 )  VFmt
C      Write( GrupFl, 18 )  Title
C      Do 1000 la = 1, LNoGap, 1
C         l = NoGaps(la)
C         Do 990 i = 1, GrpKnt-1, 1
C            TmpEnt(i) = NInt( 100.0 * ( Cross(l,1,i)+Cross(l,2,i) ))
C            TmpCh(i) = BestSq(l,1,i)
C  990       continue
C         Write( *, 19 )  VFmt
C         Write( *, * )  '  GrpKnt - 1 =',  GrpKnt-1
C         Write( GrupFl, VFmt ) l, NInt( 100.0 * Entropy(l)),
C     *               BestSqT(l),
C     1               (( NInt( 100.0 * ( Cross(l,1,i)+Cross(l,2,i))),
C     1               BestSq(l,1,i)), i = 1, GrpKnt-1, 1 )
C 1000    continue
C      Close( Unit = GrupFl, Status = 'keep' )
C
      Return
      End
C
      Subroutine WrPlot( GrupFl, GrpKnt, NotID, LGpNam, GrpNam, Title,
     *                   LNoGap, NoGaps, BestSqT, Entropy, BestSq,
     1                   Cross, AlOthr, NSP, MxGrpP, MaxSqP )
C
C ***  Write the family and group entropies, as well as consensus
C ***   sequences to a file formatted for Psi-Plot
C
      Implicit None
C
C ***   Passed Variable Declarations
C
      Integer        GrupFl, NSP, MxGrpP, MaxSqP, NotID, LNoGap, GrpKnt
      Integer        LGpNam( MxGrpP ), NoGaps( MaxSqP )
      Real*4         Entropy( MaxSqP ), Cross( MaxSqP, 2, MxGrpP ),
     *               AlOthr( MaxSqP, 2, MxGrpP )
      Character*1    BestSq( MaxSqP, 2, MxGrpP ), BestSqT( MaxSqP )
      Character*30   GrpNam( MxGrpP )
      Character*70   Title
C
C ***   Local Variable Declarations
C
      Integer        l, la, DoGrup
      Character*34   CVFile
C
    1 Format( '//')
    2 Format( '//', A70 )
    3 Format( '// Group: ', I4, ' ... ', A30 )
   17 Format( '//  Alignment', 9X, 'Group: ', A30 )
C   18 Format( 'Index  SeqAln  Entropy SeqGrp  GroupEntropy  PartGroup',
C     *        ' SeqG  outGroup    OtherDis  Other1   Other2' )
   19 Format( 'Index  SeqAln  Entropy SeqGrp  GroupEntropy  PartGroup',
     *        ' SeqNotGroup' )
C   20 Format( I5, ',  ', A1, ', ', F8.3,', ', A1, ', ', F8.3, ',',
C     *        F8.3, ', ', A1, ', ',F8.3,',', F8.3, ',', F8.3,',', F8.3)
   21 Format( I5, ',  ', A1, ', ', F8.3,', ', A1, ', ', F8.3, ',',
     *        F8.3, ', ', A1 )
C   22 Format( '//Field  1:  Alignment position index',
C     * /'//Field  2:  High entropy residue for the entire alignment',
C     1 /'//Field  3:  Entropy, at position, for the entire alignment',
C     2 /'//Field  4:  High entropy residue for the Group',
C     3 /'//Field  5:  Total Entropy for the Group (distance)',
C     4 /'//Field  6:  Partial Entropy for the Group {p(i)/q(i)}',
C     5 /'//Field  7:  High entropy residue for the .not.Group',
C     6 /'//Field  8:  Partial Entropy for the .not.Group {q(i)/p(i)}',
C     7 /'//Field  9:  Total group Entropy summed over All Other Groups',
C     8 /'//Field 10:  Partial Entropy summed over All Other Groups',
C     9               '{p(i)/q(i)}',
C     A /'//Field 11:  Partial Entropy summed over All Other Groups',
C     B               '{q(i)/p(i)}' )
   23 Format( '//Field  1:  Alignment position index',
     * /'//Field  2:  High entropy residue for the entire alignment',
     1 /'//Field  3:  Entropy, at position, for the entire alignment',
     2 /'//Field  4:  High entropy residue for the Group',
     3 /'//Field  5:  Total Entropy for the Group (distance)',
     4 /'//Field  6:  Partial Entropy for the Group {p(i)/q(i)}',
     5 /'//Field  7:  High entropy residue for the .not.Group' )
   24 Format( ' ',  4F10.3 )
   25 Format(  ' Family Entropy',
     1        /' Group Entropy',
     3        /' Partial Group Entropy',
     4        /' Partial Out Group Entropy' )
   26 Format(' ', I4, A1,'-',A1 )
C
C
      Do 300 DoGrup = 1, GrpKnt, 1
         If( DoGrup .eq. NotID )    GoTo 300
         CVFile = '                                  '
         CVFile(1:LGpNam(DoGrup)) = GrpNam(DoGrup)(1:LGpNam(DoGrup))
         CVFile(LGpNam(DoGrup)+1:LGpNam(DoGrup)+4) = '.plt'
         Open( Unit = GrupFl, File = CVFile, Status = 'new',
     *         Access = 'Sequential', form = 'Formatted')
c     *         CarriageControl = 'List' )
         Write( GrupFl, 1 )
         Write( GrupFl, 2 )   Title
         Write( GrupFl, 1 )
         Write( GrupFl, 3 )   DoGrup, GrpNam(DoGrup)
         Write( GrupFl, 1 )
C         Write( GrupFl, 22 )
         Write( GrupFl, 23 )
         Write( GrupFl, 1 )
         Write( GrupFl, 1 )
         Write( GrupFl, 17 )  GrpNam(DoGrup)
C         Write( GrupFl, 18 )
         Write( GrupFl, 19 )
         Do 200 la = 1, LNoGap, 1
            l = NoGaps(la)
C            Write( GrupFl, 20 )  l, BestSqT(l), Entropy(l),
C     *             BestSq(l,1,DoGrup),
C     1             Cross(l,1,DoGrup) + Cross(l,2,DoGrup),
C     2             Cross(l,1,DoGrup), BestSq(l,2,DoGrup),
C     3             Cross(l,2,DoGrup),
C     4             AlOthr(l,1,DoGrup) + AlOthr(l,2,DoGrup),
C     5             AlOthr(l,1,DoGrup), AlOthr(l,2,DoGrup)
            Write( GrupFl, 21 )  l, BestSqT(l), Entropy(l),
     *             BestSq(l,1,DoGrup),
     1             Cross(l,1,DoGrup) + Cross(l,2,DoGrup),
     2             Cross(l,1,DoGrup), BestSq(l,2,DoGrup)
  200       continue
         Close( Unit = GrupFl, Status = 'keep' )
  300    continue
C
C ***   Write files for XGobi plotting program.  A data (.dat) file
C ***   containing the family entropy and the group entropy as well
C ***   as the partial group entropy and the partial .not.group
C ***   entropy.  A columns definition file (.col) that labels the
C ***   four columns of data in the data file.  A rows definition (.row)
C ***   file that provides a label for each data point.  Each label
C ***   will contain the alignment index, the group residue and the
C ***   family residue.
C
      Do 800 DoGrup = 1, GrpKnt, 1
         If( DoGrup .eq. NotID )    GoTo 800
         CVFile = '                                  '
         CVFile(1:LGpNam(DoGrup)) = GrpNam(DoGrup)(1:LGpNam(DoGrup))
         CVFile(LGpNam(DoGrup)+1:LGpNam(DoGrup)+4) = '.dat'
         Open( Unit = GrupFl, File = CVFile, Status = 'new',
     *         Access = 'Sequential', form = 'Formatted')
c     *         CarriageControl = 'List' )
         Do 500 la = 1, LNoGap, 1
            l = NoGaps(la)
            Write( GrupFl, 24 )  Entropy(l),
     *                           Cross(l,1,DoGrup) + Cross(l,2,DoGrup),
     1                           Cross(l,1,DoGrup), Cross(l,2,DoGrup)
  500       continue
         Close( Unit = GrupFl, Status = 'keep' )
         CVFile = '                                  '
         CVFile(1:LGpNam(DoGrup)) = GrpNam(DoGrup)(1:LGpNam(DoGrup))
         CVFile(LGpNam(DoGrup)+1:LGpNam(DoGrup)+4) = '.row'
         Open( Unit = GrupFl, File = CVFile, Status = 'new',
     *         Access = 'Sequential', form = 'Formatted')
C     *         CarriageControl = 'List' )
         Do 600 la = 1, LNoGap, 1
            l = NoGaps(la)
            Write( GrupFl, 26 )  l, BestSq(l,1,DoGrup), BestSqT(l)
  600       continue
         Close( Unit = GrupFl, Status = 'keep' )
         CVFile = '                                  '
         CVFile(1:LGpNam(DoGrup)) = GrpNam(DoGrup)(1:LGpNam(DoGrup))
         CVFile(LGpNam(DoGrup)+1:LGpNam(DoGrup)+4) = '.col'
         Open( Unit = GrupFl, File = CVFile, Status = 'new',
     *         Access = 'Sequential', form = 'Formatted')
c     *         CarriageControl = 'List' )
         Write( GrupFl, 25 )
         Close( Unit = GrupFl, Status = 'keep' )
  800    continue
C
      Return
      End
C
      Subroutine WAsnGp( GrupFl, Title, GrpNam, GrpKnt, NGroup, InvGrp,
     *                   SqName, Fscore, GrpAln, CVScor, GrupId, NS,
     1                   NSP, MxGrpP )
C
C ***   Write Assignment to Groups information to a file
C
C ***   Passed Variable declarations
C
      Integer       NSP, MxGrpP, NS, GrupFl, GrpKnt
      Integer       Ngroup( MxGrpP ), InvGrp( NSP, MxGrpP ), GrupId(NSP)
      Real*4        CVScor( NSP, 2 ), GrpAln( NSP, 0:MxGrpP ),
     *              Fscore( NSP )
      Character*15  SqName( NSP )
      Character*30  GrpNam( MxGrpP )
      Character*70  Title
C
C ***   Local varaible declarations
C
      Integer      DoGrup, Nblock, ib, ie, n, sn
C
    2 Format(/' ', A70, / )
    3 Format(' Group: ', I4, ' ... ', A30 )
   21 Format(/'  Sequence Name  FamilyScore  ', 8( 'Group', I3:,3X))
C   22 Format(' ', 29X, 8( 'Group', I3:,3X))
   23 Format( A15, 2X, F10.3, 3X, 8F10.3: )
   24 Format( ' ', 29X, 8F10.3: )
   25 Format(//' Sequence Name  FamilyScores  GroupNumber',
     *         '  GroupScores',/)
   26 Format( A15, 1X, 3F10.3, 1X, I3, 1X, 8F10.3: )
   27 Format( ' ', 51X, 8F10.3: )
C
C ***   Write the report file that assigns previously unclassified
C ***   sequences to specific groups.  Open the file and list the
C ***   user defined groups to which sequences can be assigned.
C
      Open( Unit = GrupFl, File = 'AssignSequences.cen',
     *      Access = 'Sequential', form = 'Formatted',
     *      Status = 'new' )
c     *      CarriageControl = 'List', Status = 'new' )
      Write( GrupFl, 2 )  Title
      Write( GrupFl, 3 )  ( DoGrup, GrpNam( DoGrup ), DoGrup = 1,
     *                      GrpKnt-1, 1 )
      NBlock = ( GrpKnt + 6 ) / 8
      ib = 1
      ie = 8
      If( ie .gt. GrpKnt-1 )    ie = GrpKnt-1
      Write( GrupFl, 21 )  ( DoGrup, DoGrup = ib, ie, 1 )
      If( NBlock .gt. 1 )    Then
         Do 1300 nb = 2, NBlock, 1
            ib = ib + 8
            ie = ib + 7
            If( ie .gt. GrpKnt-1 )    ie = GrpKnt-1
            Write( GrupFl, 21 )  ( DoGrup, DoGrup = ib, ie, 1 )
 1300       continue
      Else
      EndIf
C
C ***   1500 loop over all sequences in the highest numbered formal group.
C ***   This group contains the sequences designated by the user to be
C ***   assigned to one of the user defined groups of sequences.
C
      Do 1500 n = 1, NGroup( GrpKnt ), 1
         sn = InvGrp(n,GrpKnt)
         ib = 1
         ie = 8
         If( ie .gt. GrpKnt-1 )    ie = GrpKnt-1
         Write( GrupFl, 23 ) SqName(sn), FScore(sn),
     *                       ( GrpAln(sn,DoGrup), DoGrup = ib, ie, 1 )
         If( NBlock .gt. 1 )    Then
            Do 1400 nb = 2, NBlock, 1
               ib = ib + 8
               ie = ib + 7
               If( ie .gt. GrpKnt-1 )    ie = GrpKnt-1
               Write( GrupFl, 24 ) ( GrpAln(sn,DoGrup),DoGrup=ib,ie,1)
 1400          continue
          Else
          EndIf
 1500    continue
C
C ***   Report alignment scores for all sequences with all groups
C
      Write( GrupFl, 25 )
      NBlock = ( GrpKnt + 6 ) / 8
      Do 1700 n = 1, NS, 1
         ib = 1
         ie = 8
         If( ie .gt. GrpKnt-1 )    ie = GrpKnt-1
         Write( GrupFl, 26 )  SqName(n), FScore(n), CVScor(n,1),
     *                        CVScor(n,2), GrupId(n),
     1                      ( GrpAln(n,DoGrup), DoGrup = ib, ie, 1 )
         If( NBlock .gt. 1 )    Then
            Do 1600 nb = 2, NBlock, 1
               ib = ib + 8
               ie = ib + 7
               If( ie .gt. GrpKnt-1 )    ie = GrpKnt-1
               Write( GrupFl, 27 ) ( GrpAln(n,DoGrup),DoGrup=ib,ie,1)
 1600          continue
          Else
          EndIf
 1700    continue
      Close( Unit = GrupFl, Status = 'keep' )
C
      Return
      End
C
C
      Subroutine WrPSSM( Termnl, GrupFl, ProfileT, Profil, MaxSqP,
     *                   MxGrpP, AlfMap, MapAA, MapNuc, NoGaps, MemeAA,
     1                   MemNuc, Alfbet, GrpNam, Title, Protin, PFile,
     2                   ASize, LNoGap, thresh, GrpKnt, PFLen, AlnLen,
     3                   MxSymP, Comp )
C
C ***   Write Possition Specific Scoring Matrices (profiles) for the
C ***   family and each group to a file
C
      Implicit None
C
C ***   Passed variable declarations
C
      Integer         Termnl, GrupFl, MaxSqP, MxGrpP, MxSymP
      Integer         AlfMap(20), MapAA( 20 ), MapNuc(4), PFLen,
     *                NoGaps(MaxSqP), LNoGap, ASize, GrpKnt, AlnLen,
     1                Comp(0:MxSymP)
      Real*4          Profil( 25, MaxSqP, 2, MxGrpP ),
     *                ProfileT(23, MaxSqP ), thresh( MxGrpP, 2 )
      Character*1     MemeAA(20), MemNuc(4), Alfbet(20)
      Character*24    PFile
      Character*30    GrpNam( MxGrpP )
      Character*70    Title
      Logical         Protin
C
C ***   Local variable declarations
C
      Integer         i, l, la, DoGrup, TotCh
      Character*34    CVFile
C
    2 Format(/' ', A70, / )
   28 Format(//' Meme Style log-odds Position Specific Scoring',
     *         ' Matrices.',
     1       //' Residue Order conforms to MEME output and is',
     2         ' alphabetical.',
     3        /'   ', 20( 1X, A1: ) )
   29 Format(//' The first matrix is for the entire family.', // )
   30 Format( 20( F8.3: ) )
   31 Format(//'log-odds matrix: alength=', I3,' w=', I3,' n=', I6,
     *         ' bayes=', F11.4 )
   32 Format(/' Position Specific Scoring Matrix for recognizing',
     *        ' members of Group:',
     1       /'      Group Number =', I3, ',   Name: ', A30, /)
   33 Format(/' Position Specific Scoring Matrix for recognizing',
     *        ' not-members of Group:',
     1       /'      Group Number =', I3, ',   Name: ', A30, /)
C
C ***   Write the profiles into a file
C
      CVFile = '                                  '
      CVFile(1:PfLen) = PFile(1:PfLen)
      CVFile(PfLen+1:PfLen+4) = '.prf'
C      Write( Termnl, 3 )  PfLen, PFile, CVFile                            DeBug
C    3 Format(/' PfLen =', I6, ',   PFile = ', A24,                        DeBug
C     *       /' and CVFile = ', A34 )                                     DeBug
      Open( Unit = GrupFl, File = CVFile, Status = 'new',
     *      Access = 'Sequential', form = 'Formatted',
     *      RECL = 164 )
C     *      CarriageControl = 'List', RECL = 164 )
C
C      Write( Termnl, 10 )                                                 DeBug
C   10 Format(/' PSSM file opened.' )                                      DeBug
      TotCh = 0
      Do 1800 l = 1, ASize, 1
         TotCh = TotCh + Comp(l)
 1800    continue
C      Write( Termnl, 11 )  TotCh                                          DeBug
C   11 Format(/' Total Characters =', I15 )                                DeBug
      If( Protin )    Then
         Do 1820 l = 1, ASize, 1
            Alfbet(l) = MemeAA(l)
            AlfMap(l) = MapAA(l)
 1820       continue
      Else
         Do 1840 l = 1, ASize, 1
            Alfbet(l) = MemNuc(l)
            AlfMap(l) = MapNuc(l)
 1840       continue
      EndIf
C      Write( Termnl, 2 )  Title
      Write( GrupFl, 2 )  Title
      Write( GrupFl, 28 ) ( Alfbet(l), l = 1, ASize, 1 )
      Write( GrupFl, 29 )
      Write( GrupFl, 31 )  ASize, AlnLen, TotCh, 99.9999
C      Write( Termnl, 28 ) ( Alfbet(l), l = 1, ASize, 1 )
      Do 1900 la = 1, LNoGap, 1
         l = NoGaps(la)
CAJR added next 3 for debug 2020:
C         Write( 6, *) l
C         Write( 6, *)  ( AlfMap(i), i =1,ASize,1)
C         Write( 6, *)  ( ProfileT( AlfMap(i), l ), i =1,ASize,1)
         Write( GrupFl, 30)  ( ProfileT( AlfMap(i), l ), i =1,ASize,1)
 1900    continue
C      Write( Termnl, 12 )                                                 DeBug
C   12 Format(/' Total Alignment PSSM written.' )                          DeBug
      Do 2200 DoGrup = 1, GrpKnt-1, 1
         Write( GrupFl, 2 )  Title
         Write( GrupFl, 32 )  DoGrup, GrpNam( DoGrup )
C         Write( Termnl, 32 )  DoGrup, GrpNam( DoGrup )                    DeBug
         Write( GrupFl, 31 )  ASize, AlnLen, TotCh, thresh(DoGrup,1)
         Do 2000 la = 1, LNoGap, 1
            l = NoGaps(la)
            Write(GrupFl,30) (Profil(AlfMap(i),l,1,DoGrup),i=1,ASize,1)
 2000       continue
         Write( GrupFl, 33 )  DoGrup, GrpNam( DoGrup )
         Write( GrupFl, 31 )  ASize, AlnLen, TotCh, thresh(DoGrup,2)
C         Write( Termnl, 31 )  ASize, AlnLen, TotCh, thresh(DoGrup,2)      DeBug
         Do 2100 la = 1, LNoGap, 1
            l = NoGaps(la)
            Write(GrupFl,30) (Profil(AlfMap(i),l,2,DoGrup),i=1,ASize,1)
 2100       continue
 2200    continue
      Close( Unit = GrupFl, Status = 'keep' )
C
      Return
      End
C
C
      Subroutine EWhole( NSP, TotSqP, MaxSqP, n20, NS, Seq, LsCum,
     *                   GrupID, ProfileT, Qij, QSum, Entropy, ASize,
     1                   AlnLen, NoGaps, LNoGap, Protin, MPC, SwProt,
     2                   BestSqT, MaxLet, SqPrnt, alnEnt, AlnScr,
     3                   FScore, DeBugF, NotID, Consen )
C
C ***   Compute the entropy for each position in the family
C
      Implicit None
C
C ***   Passed Variables - Subroutine EWhole (Entropy calculations)
C
      Integer       NSP, TotSqP, MaxSqP, n20, NS, MaxLet, ASize,
     *              AlnLen, DeBugF, LNoGap, NotID
      Integer       LsCum( NSP ), GrupId( NSP ), Seq( TotSqP ),
     *              NoGaps( MaxSqP )
      Real*4        ProfileT( 23, MaxSqP ), Entropy( MaxSqP ), MPC,
     *              Qij( n20, n20 ), QSum( n20 ), SwProt( n20 ),
     1              AlnEnt, AlnScr, FScore( NSP )
      Character*1   BestSqT( MaxSqP ), Consen( MaxSqP ), sqprnt(25)
      Logical       Protin
C
C ***   Local variables
C
      Integer       i, la, l, n, m, bestt, aatype, high
      Real*4        Tknt, ttot, t(25), tfake(25), MaxT, HiKnt
      Real*8        base2, odds
C
C
      base2 = 1.0D 00 / dlog( 2.0D 00 )
      AlnEnt = 0.0
      AlnScr = 0.0
      Do 100 i = 1, NSP, 1
         FScore(i) = 0.0
  100    continue
      Do 1000 la = 1, LNoGap, 1
         l = NoGaps(la)
         Tknt = 0.0
         ttot = 0.0
         Entropy(l) = 0.0
         Do 160 i = 1, 25, 1
            t(i) = 0.0
            tfake(i) = 0.0
  160       continue
C
C ***   Get the raw counts of sequence residues for all sequences and
C ***   at the current position
C
         Do 200 n = 1, NS, 1
            If( GrupId(n) .ne. 0 .and. GrupId(n) .ne. NotID )    Then
               aatype = Seq(LsCum(n)+l)
               t( aatype ) = t( aatype ) + 1.0
               If( aatype .ge. 1  .and.  aatype .le. MaxLet )
     *                Ttot = Ttot + 1.0
            Else
            EndIf
  200       continue
C
C ***   If the sequences are protein adjust the counts for the ambiguous
C ***   amino acids Asx = (Asn,Asp) and Glx = (Gln,Glu).
C
         If( Protin )    Then
            t(3) = t(3) + ( t(21) / 2.0 )
            t(4) = t(4) + ( t(21) / 2.0 )
            t(21) = 0.0
            t(6) = t(6) + ( t(22) / 2.0 )
            t(7) = t(7) + ( t(22) / 2.0 )
            t(22) = 0.0
         Else
         EndIf
C    9 Format(/'  Total (Ttot) number of amino acids in column L =',
C     *           F6.1)
C   10 Format('  Raw values for t(n), n = 1 -> 20;  L = ', I3 )
C   11 Format( 5X,  10F7.1 )
C      Write( DeBugF, 9 )   Ttot
C      Write( DeBugF, 10 )  l
C      Write( DeBugF, 11 ) ( t(n), n =  1, 10, 1 )
C      Write( DeBugF, 11 ) ( t(n), n = 11, 20, 1 )
C
C ***   Count the different kinds of sequence residue present at this
C ***   position and multiply it by MPC, the pseudocount multiplier.
C ***   Set the consensus sequence residue to the residue type with the
C ***   highest count.
C
         HiKnt = 0.0
         high = 0
         Do 250 n = 1, ASize, 1
            If( t(n) .gt. HiKnt )    Then
               HiKnt = t(n)
               high = n
            Else
            EndIf
            If( t(n) .ge. 0.0001 )    Tknt = Tknt + MPC
  250       continue
         Consen(l) = SqPrnt(high)
C
C ***   Compute the Henikoff style pseudocounts to be added to each
C ***   kind of sequence residue at this position.
C
         Do 310 n = 1, ASize, 1
            Do 300 m = 1, ASize, 1
C2020 Debug Statement
C              write(6,*) 'TFAKE:',tfake(n),t(m),Ttot,Qij(m,n),Qsum(m)
               tfake(n) = tfake(n) + ((t(m)/Ttot) * (Qij(m,n)/QSum(m)))
  300          continue
            tfake(n) = tfake(n) * Tknt
  310       continue
C    8 Format(/'  Total pseudocounts for this position =', F6.1 )
C   12 Format('  Raw tfake(n), n = 1 -> 20;  L = ', I3 )
C   13 Format( 5X,  10F7.4 )
C      Write( DeBugF, 8 )   Tknt
C      Write( DeBugF, 12 )  l
C      Write( DeBugF, 13 ) ( tfake(n), n =  1, 10, 1 )
C      Write( DeBugF, 13 ) ( tfake(n), n = 11, 20, 1 )
C
C ***   Compute the corrected counts as the weighted average of the
C ***   actual, observed counts and the computed pseudocounts
C
         Do 330 n = 1, ASize, 1
C2020 Debug Print
C            write(6,*) Ttot, Tknt, t(n), tfake(n)
            t(n) = (( Ttot / ( Ttot + Tknt )) * ( t(n) / Ttot ))
     *           + (( Tknt / ( Ttot + Tknt )) * ( tfake(n) / Tknt ))
  330       continue
C   14 Format(/'  Supplemented t(n), n = 1 -> 20;  L = ', I3 )
C      Write( DeBugF, 14 )  l
C      Write( DeBugF, 13 ) ( t(n), n =  1, 10, 1 )
C      Write( DeBugF, 13 ) ( t(n), n = 11, 20, 1 )
C
C ***   Compute the entropy for all of the aligned sequences.  Then
C ***   compute the profile or Position Specific log-odds Scoring Matrix
C ***   in bits.  Also find the sequence that has the best score for
C ***   the profile.
C
         MaxT = -1000000.0
         bestt = 0
         Do 350 n = 1, ASize, 1
            odds = t(n) / SwProt(n)
            If( odds .le. 0.000001 )   Write( *, * )  ' L=', l,           DeBug
     *           ',   n=', n, ',   t(n) =', t(n)                          DeBug
C 2020 print out for debug:
C            write(6,*) base2, t(n), SwProt(n), odds, 
C     +            dlog(odds), base2 * dlog(odds)
            ProfileT(n,l) = base2 * dlog( odds )
            Entropy(l) = Entropy(l) + ( t(n) * ProfileT(n,l) )
            If( ProfileT(n,l) .gt. MaxT )    Then
               MaxT = ProfileT(n,l)
               bestt = n
            Else
            EndIf
  350       continue
         BestSqT(l) = sqprnt(bestt)
         AlnScr = AlnScr + ProfileT(bestt,l)
         AlnEnt = AlnEnt + Entropy(l)
C   15 Format(/' Family Profile(n), n = 1 -> 20; for L =', I4,            DeBug
C     *        ',   Best = ', A1 )                                        DeBug
C      Write( DeBugF, 15 )  l, BestSqT(l)                                 DeBug
C      Write( DeBugF, 13 ) ( ProfileT(n,l), n =  1, 10, 1 )               DeBug
C      Write( DeBugF, 13 ) ( ProfileT(n,l), n = 11, 20, 1 )               DeBug
C
C ***   Compute B=Asx (21), Z=Glx (22), and X (23) averages and scores.
C
         If( Protin )    Then
            ProfileT(21,l) = (( ProfileT(3,l) * t(3) )
     *                     +  ( ProfileT(4,l) * t(4) )) / (t(3)+t(4))
            ProfileT(22,l) = (( ProfileT(6,l) * t(6) )
     *                     +  ( ProfileT(7,l) * t(7) )) / (t(6)+t(7))
         Else
         EndIf
         ProfileT(MaxLet,l) = 0.0
         MaxT = 0.0
         Do 400 n = 1, ASize, 1
            ProfileT(MaxLet,l) = ProfileT(MaxLet,l)
     *                         + (ProfileT(n,l)* t(n))
            MaxT = MaxT + t(n)
  400       continue
         ProfileT(MaxLet,l) = ProfileT(MaxLet,l) / MaxT
C
C ***   Score all sequences in the alignment against the family profile
C
         Do 500 n = 1, NS, 1
            FScore(n) = FScore(n) + ProfileT(Seq(LsCum(n)+l),l)
  500       continue
C
 1000    continue
C
      Return
      End
C
C
      Subroutine GrpEnt( NSP, TotSqP, MaxSqP, n20, NS, MxSymP, Seq,
     *                   LsCum, GrupID, Qij, QSum, Profil, Cross,
     1                   RawKnt, BackF, AlnLen, NoGaps, LNoGap, ASize,
     2                   Protin, MPC, BestSq, GCons, Gap, MaxLet,
     3                   SqPrnt, MxGrpP, NotID, DoGrup, MGpScr, GrpAln,
     4                   DeBugF, TotCE )
C
C ***   Computes cross-entropies for each group.  The .not.group consist
C ***   of all of the sequence that have been assigned to another
C ***   analysis group.  If a sequence is in the "ToBeClassified" group
C ***   or in the group of sequences to ignore (Group Number = zero)
C ***   it is not included in any analysis group or in the .not.group
C ***   for any analysis group.  The passed variable DoGrup identifies
C ***   the group number (array GrupId) of the current analysis group.
C ***   The passed variable NotID identifies the group number of the
C ***   sequence that are to be classified later.
C
      Implicit None
C
C ***   Passed Variables - Subroutine GrpEnt
C
      Integer       NSP, TotSqP, MaxSqP, n20, NS, MxSymP, NotID, DeBugF,
     *              AlnLen, ASize, MxGrpP, DoGrup, MaxLet, Gap, LNoGap
      Integer       LsCum( NSP ), GrupId( NSP ), Seq( TotSqP ),
     *              NoGaps( MaxSqP )
      Real*4        Profil(25, MaxSqP, 2, MxGrpP ), QSum( n20 ),
     *              Cross( MaxSqP, 2, MxGrpP ), MPC, Qij( n20, n20 ),
     1              RawKnt(25,MaxSqP), BackF(25,MaxSqP),
     2              GrpAln( NSP, 0:MxGrpP ), TotCE( MxGrpP, 2 ),
     3              MGpScr( MxGrpP )
      Character*1   BestSq( MaxSqP, 2, MxGrpP ), sqprnt(25),
     *              GCons( MaxSqP, 2, MxGrpP )
      Character*30  GrpNam( MxGrpP )
      Logical       Protin
C
C ***   Local variables
C
      Integer       i, la, l, n, m, best1, best2, aatype, pmi, qmi
      Real*4        Pknt, Qknt, ptot, qtot, p(25), q(25), pfake(25),
     *              qfake(25), Max1, Max2, pgaps, qgaps, pmax, qmax
      Real*8        base2, odds
C
C
      base2 = 1.0D 00 / dlog( 2.0D 00 )
      Do 100 n = 1, NS, 1
         GrpAln( n, DoGrup ) = 0.0
  100    continue
      TotCE( DoGrup, 1 ) = 0.0
      TotCE( DoGrup, 2 ) = 0.0
      MGpScr( DoGrup ) = 0.0
C
C ***   Loop 1000 is a loop over the number of alignment positions that
C ***   do not have gaps in any of the sequences that belong either to
C ***   analysis groups or to the "ToBeClassified Group".  The array
C ***   NoGaps is a list of those positions created in Subroutine RmGaps.
C
      Do 1000 la = 1, LNoGap, 1
         l = NoGaps(la)
C
C      Do 1000 l = 1, AlnLen, 1
         Pknt = 0.0
         Qknt = 0.0
         ptot = 0.0
         qtot = 0.0
         pgaps = 0.0
         qgaps = 0.0
         pmax = 0.0
         qmax = 0.0
         pmi = 0
         qmi = 0
         Cross( l, 1, DoGrup ) = 0.0
         Cross( l, 2, DoGrup ) = 0.0
         Do 150 i = 1, 25, 1
            p(i) = 0.0
            q(i) = 0.0
            pfake(i) = 0.0
            qfake(i) = 0.0
  150       continue
C
C ***   Get the raw counts of sequence residues for all sequences in
C ***   the analysis group and the .not.group at the current position.
C
         Do 200 n = 1, NS, 1
            If( GrupId(n) .eq. DoGrup )    Then
               aatype = Seq(LsCum(n)+l)
               p( aatype ) = p( aatype ) + 1.0
               If( aatype .ge. 1  .and.  aatype .le. MaxLet )    Then
                  Ptot = Ptot + 1.0
               ElseIf( aatype .eq. Gap )    Then
                  pgaps = pgaps + 1.0
               Else
               EndIf
            ElseIf( GrupId(n) .gt. 0   .and.
     *              GrupId(n) .ne. NotID )    Then
               aatype = Seq(LsCum(n)+l)
               q( aatype ) = q( aatype ) + 1.0
               If( aatype .ge. 1  .and.  aatype .le. MaxLet )    Then
                  Qtot = Qtot + 1.0
               ElseIf( aatype .eq. Gap )    Then
                  qgaps = qgaps + 1.0
               Else
               EndIf
            Else
            EndIf
  200       continue
C
C ***   Accumulate the raw counts for the foreground group in an array
C ***   over all positions in the alignment to use for cross-validation
C ***   of the members of the foreground group,
C
         Do 210 i = 1, Gap, 1
            RawKnt(i,l) = p(i)
  210       continue
C
C ***   If the sequences are protein adjust the counts for the ambiguous
C ***   amino acids Asx = (Asn,Asp) and Glx = (Gln,Glu).
C
         If( Protin )    Then
            p(3) = p(3) + ( p(21) / 2.0 )
            p(4) = p(4) + ( p(21) / 2.0 )
            p(21) = 0.0
            p(6) = p(6) + ( p(22) / 2.0 )
            p(7) = p(7) + ( p(22) / 2.0 )
            p(22) = 0.0
            q(3) = q(3) + ( q(21) / 2.0 )
            q(4) = q(4) + ( q(21) / 2.0 )
            q(21) = 0.0
            q(6) = q(6) + ( q(22) / 2.0 )
            q(7) = q(7) + ( q(22) / 2.0 )
            q(22) = 0.0
         Else
         EndIf
C
C ***   Find the most common sequence residue at this position
C
         Do 220 i = 1, ASize, 1
            If( p(i) .gt. pmax )    Then
               pmax = p(i)
               pmi = i
            Else
            EndIf
            If( q(i) .gt. qmax )    Then
               qmax = q(i)
               qmi = i
            Else
            EndIf
  220       continue
         If( pgaps .lt. Ptot )    Then
            GCons(l,1,DoGrup) = sqprnt(pmi)
         Else
            GCons(l,1,DoGrup) = sqprnt(Gap)
         EndIf
         If( qgaps .lt. qtot )    Then
            GCons(l,2,DoGrup) = sqprnt(qmi)
         Else
            GCons(l,2,DoGrup) = sqprnt(Gap)
         EndIf
C    7 Format(/' Indirect Index, L=', I4, ',   LA=', I4,                  DeBug
C     *        ',   NoGaps(LA)=', I4 )                                    DeBug
C    9 Format(/'  Ptot = number of Group amino acids in column L =',      DeBug
C     *           F6.1)                                                   DeBug
C   10 Format('  Raw values for p(n), n = 1 -> 20;  L = ', I3 )           DeBug
C   29 Format(/'  Qtot = number .not.Group amino acids in column L =',    DeBug
C     *           F6.1)                                                   DeBug
C   30 Format('  Raw values for q(n), n = 1 -> 20;  L = ', I3 )           DeBug
C   11 Format( 5X,  10F7.1 )                                              DeBug
C      Write( DeBugF, 7 )   l, la, NoGaps(la)                             DeBug
C      Write( DeBugF, 9 )   Ptot                                          DeBug
C      Write( DeBugF, 10 )  l                                             DeBug
C      Write( DeBugF, 11 ) ( p(n), n =  1, 10, 1 )                        DeBug
C      Write( DeBugF, 11 ) ( p(n), n = 11, 20, 1 )                        DeBug
C      Write( DeBugF, 29 )   Qtot                                         DeBug
C      Write( DeBugF, 30 )  l                                             DeBug
C      Write( DeBugF, 11 ) ( q(n), n =  1, 10, 1 )                        DeBug
C      Write( DeBugF, 11 ) ( q(n), n = 11, 20, 1 )                        DeBug
C
C ***   Count the different kinds of sequence residue present in each
C ***   group at this position and multiply it by the pseudocount
C ***   multiplier, MPC.
C
         Do 350 n = 1, ASize, 1
            If( p(n) .ge. 0.0001 )    Pknt = Pknt + MPC
            If( q(n) .ge. 0.0001 )    Qknt = Qknt + MPC
  350       continue
C
C ***   Compute the Henikoff style pseudocounts to be added to each
C ***   kind of sequence residue at this position.
C
         Do 410 n = 1, ASize, 1
            Do 400 m = 1, ASize, 1
               pfake(n) = pfake(n) + ((p(m)/Ptot) * (Qij(m,n)/QSum(m)))
               qfake(n) = qfake(n) + ((q(m)/Qtot) * (Qij(m,n)/QSum(m)))
  400          continue
            pfake(n) = pfake(n) * Pknt
            qfake(n) = qfake(n) * Qknt
  410       continue
C    8 Format(/'  Total, Pknt, pseudocounts for this position =', F6.1)   DeBug
C   12 Format('  Raw pfake(n), n = 1 -> 20;  L = ', I3 )                  DeBug
C   28 Format(/'  Total, Qknt, pseudocounts for this position =', F6.1)   DeBug
C   32 Format('  Raw qfake(n), n = 1 -> 20;  L = ', I3 )                  DeBug
C   13 Format( 5X,  10F7.4 )                                              DeBug
C      Write( DeBugF, 8 )   Pknt                                          DeBug
C      Write( DeBugF, 12 )  l                                             DeBug
C      Write( DeBugF, 13 ) ( pfake(n), n =  1, 10, 1 )                    DeBug
C      Write( DeBugF, 13 ) ( pfake(n), n = 11, 20, 1 )                    DeBug
C      Write( DeBugF, 28 )  Qknt                                          DeBug
C      Write( DeBugF, 32 )  l                                             DeBug
C      Write( DeBugF, 13 ) ( qfake(n), n =  1, 10, 1 )                    DeBug
C      Write( DeBugF, 13 ) ( qfake(n), n = 11, 20, 1 )                    DeBug
C
C ***   Compute the corrected counts as the weighted average of the
C ***   actual, observed counts and the computed pseudocounts
C
         Do 430 n = 1, ASize, 1
            p(n) = (( Ptot / ( Ptot + Pknt )) * ( p(n) / Ptot ))
     *           + (( Pknt / ( Ptot + Pknt )) * ( pfake(n) / Pknt ))
            q(n) = (( Qtot / ( Qtot + Qknt )) * ( q(n) / Qtot ))
     *           + (( Qknt / ( Qtot + Qknt )) * ( qfake(n) / Qknt ))
  430       continue
C   14 Format(/'  Supplemented p(n), n = 1 -> 20;  L = ', I3 )            DeBug
C      Write( DeBugF, 14 )  l                                             DeBug
C      Write( DeBugF, 13 ) ( p(n), n =  1, 10, 1 )                        DeBug
C      Write( DeBugF, 13 ) ( p(n), n = 11, 20, 1 )                        DeBug
C   34 Format(/'  Supplemented q(n), n = 1 -> 20;  L = ', I3 )            DeBug
C      Write( DeBugF, 34 )  l                                             DeBug
C      Write( DeBugF, 13 ) ( q(n), n =  1, 10, 1 )                        DeBug
C      Write( DeBugF, 13 ) ( q(n), n = 11, 20, 1 )                        DeBug
C
C ***   Compute the entropy for each of
C ***   the two analysis groups.  Then compute the two profiles or
C ***   Position Specific Scoring Matrices log-odds matrices in bits.
C ***   Also find the sequences that have the best scores in for each of
C ***   the two profiles.
C
         Max1 = -1000000.0
         Max2 = -1000000.0
         best1 = 0
         best2 = 0
         Do 450 n = 1, ASize, 1
            odds = p(n) / q(n)
            Profil(n,l,1,DoGrup) = base2 * dlog( odds )
            Cross( l, 1, DoGrup ) = Cross(l,1,DoGrup)
     *                            + ( p(n) * Profil(n,l,1,DoGrup) )
            TotCE( DoGrup, 1 ) = TotCE(DoGrup,1) + Cross(l,1,DoGrup)
            If( Profil(n,l,1,DoGrup) .gt. Max1 )    Then
               Max1 = Profil(n,l,1,DoGrup)
               best1 = n
            Else
            EndIf
            odds = q(n) / p(n)
            Profil(n,l,2,DoGrup) = base2 * dlog( odds )
            Cross( l, 2, DoGrup ) = Cross(l,2,DoGrup)
     *                            + ( q(n) * Profil(n,l,2,DoGrup) )
            TotCE( DoGrup, 2 ) = TotCE(DoGrup,2) + Cross(l,2,DoGrup)
            If( Profil(n,l,2,DoGrup) .gt. Max2 )    Then
               Max2 = Profil(n,l,2,DoGrup)
               best2 = n
            Else
            EndIf
  450       continue
         MGpScr(DoGrup) = MGpScr(DoGrup) + Profil(best1,l,1,DoGrup)
         BestSq(l,1,DoGrup) = sqprnt(best1)
         BestSq(l,2,DoGrup) = sqprnt(best2)
C   15 Format(/' Family Profile(n), n = 1 -> 20; for L =', I4,            DeBug
C     *        ',   Best = ', A1 )                                        DeBug
C      Write( DeBugF, 15 )  l, BestSq(l,1,DoGrup)                         DeBug
C      Write( DeBugF, 13 ) ( Profil(n,l,1,DoGrup), n =  1, 10, 1 )        DeBug
C      Write( DeBugF, 13 ) ( Profil(n,l,1,DoGrup), n = 11, 20, 1 )        DeBug
C   35 Format(/' .not.Family Profile(n), n = 1 -> 20; for L =', I4,       DeBug
C     *        ',   Best = ', A1 )                                        DeBug
C      Write( DeBugF, 35 )  l, BestSq(l,2,DoGrup)                         DeBug
C      Write( DeBugF, 13 ) ( Profil(n,l,2,DoGrup), n =  1, 10, 1 )        DeBug
C      Write( DeBugF, 13 ) ( Profil(n,l,2,DoGrup), n = 11, 20, 1 )        DeBug
C
C ***   Compute B=Asx (21), Z=Glx (22), and X (23) averages and scores.
C
         If( Protin )    Then
            Profil(21,l,1,DoGrup) = (( Profil(3,l,1,DoGrup) * p(3) )
     *              +  ( Profil(4,l,1,DoGrup) * p(4) )) / (p(3) +p(4))
            Profil(21,l,2,DoGrup) = (( Profil(3,l,2,DoGrup) * q(3) )
     *              +  ( Profil(4,l,2,DoGrup) * q(4) )) / (q(3) +q(4))
            Profil(22,l,1,DoGrup) = (( Profil(6,l,1,DoGrup) * p(6) )
     *              +  ( Profil(7,l,1,DoGrup) * p(7) )) / (p(6) +p(7))
            Profil(22,l,2,DoGrup) = (( Profil(6,l,2,DoGrup) * q(6) )
     *              +  ( Profil(7,l,2,DoGrup) * q(7) )) / (q(6) +q(7))
         Else
         EndIf
         Profil(MaxLet,l,1,DoGrup) = 0.0
         Profil(MaxLet,l,2,DoGrup) = 0.0
         Max1 = 0.0
         Max2 = 0.0
         Do 500 n = 1, ASize, 1
            Profil(MaxLet,l,1,DoGrup) = Profil(MaxLet,l,1,DoGrup)
     *                            + ( Profil(n,l,1,DoGrup) * p(n) )
            Max1 = Max1 + p(n)
            Profil(MaxLet,l,2,DoGrup) = Profil(MaxLet,l,2,DoGrup)
     *                            + ( Profil(n,l,2,DoGrup) * q(n) )
            Max2 = Max2 + q(n)
  500       continue
         Profil(MaxLet,l,1,DoGrup) = Profil(MaxLet,l,1,DoGrup) / Max1
         Profil(MaxLet,l,2,DoGrup) = Profil(MaxLet,l,2,DoGrup) / Max2
C
C ***   Accumulate the pseudocount corrected background counts in an array
C ***   over all positions in the alignment to use for cross-validation
C ***   of the members of the foreground group.
C
         Do 550 i = 1, 25, 1
            BackF(i,l) = q(i)
  550       continue
C
C ***   Compute Alignment Scores for all of the sequences in the
C ***   complete alignment with the profiles for all of the groups.
C
         Do 600 n = 1, NS, 1
            GrpAln( n, DoGrup ) = GrpAln(n,DoGrup)
     *                          + Profil(Seq(LsCum(n)+l),l,1,DoGrup)
  600       continue
C
 1000    continue
C
      Return
      End
C
C
      Subroutine CRVal( NSP, TotSqP, MaxSqP, n20, MxSymP, SEQ, LsCum,
     *                  Group1, Qij, QSum, CVPEnt, CVPScr, CVEntr,
     1                  CVScor, RawKnt, BackF, AlnLen, NoGaps, LNoGap,
     2                  ASize, MaxLet, Protin, MPC, NGrp1, SqPrnt,
     3                  MxGrpP, thresh, DoGrup, Unknown, BackFr,
     4                  DeBugF )
C
      Implicit None
C
C ***   Passed Variables - Subroutine CRVal (Cross validation)
C
      Integer       NSP, TotSqP, MaxSqP, n20, MxSymP, DoGrup, unknown,
     *              AlnLen, ASize, MxGrpP, MaxLet, NGrp1, LNoGap,
     *              DeBugF
      Integer       LsCum( NSP ), Group1( NSP ), Seq( TotSqP ),
     *              NoGaps( MaxSqP )
      Real*4        CVPEnt( TotSqP, 2 ), CVPScr( TotSqP, 2 ),
     *              CVEntr( NSP, 2 ), CVScor( NSP, 2 ), MPC,
     1              Qij( n20, n20 ), QSum( n20 ), thresh( MxGrpP, 2 ),
     2              RawKnt(25,MaxSqP), BackF(25,MaxSqP), BackFr(MxSymP)
      Character*1   sqprnt(25)
      Logical       Protin
C
C ***   Local variables
C
      Integer       AAid, i, la, l, n, m, SeqPos, ncv, NSeqCV
      Real*4        Pknt, ptot, p(25), pfake(25), Max1, Max2,
     *              Profl1(25), Profl2(25)
      Real*8        base2, odds
C
C
      base2 = 1.0D 00 / dlog( 2.0D 00 )
      thresh(DoGrup,1) = -1000000.0
      thresh(DoGrup,2) = 1000000.0
C
C ***   If the analysis group contains less than three sequences the
C ***   cross validation is not a very useful exercise, so set the
C ***   cross validation results to zero and return for the next group.
C
      Do 200 ncv = 1, NGrp1, 1
         NSeqCV = Group1( ncv )
         CVEntr( NSeqCV, 1 ) = 0.0
         CVEntr( NSeqCV, 2 ) = 0.0
         CVScor( NSeqCV, 1 ) = 0.0
         CVScor( NSeqCV, 2 ) = 0.0
         Do 150 l = 1, AlnLen, 1
            SeqPos =  LsCum( NSeqCV ) + l
            CVPEnt( SeqPos, 1 ) = 0.0
            CVPEnt( SeqPos, 2 ) = 0.0
            CVPScr( SeqPos, 1 ) = 0.0
            CVPScr( SeqPos, 2 ) = 0.0
  150       continue
  200    continue
      If( NGrp1 .lt. 3 )    Then
         thresh(DoGrup,1) = -1.0
         thresh(DoGrup,2) = 1.0
         Return
      Else
      EndIf
C
C ***   Do the cross validation by removing one sequence at at time from
C ***   the analysis group.
C
      Do 2000 ncv = 1, NGrp1, 1
         NSeqCV = Group1( ncv )
         Do 1000 la = 1, LNoGap, 1
            l = NoGaps(la)
            SeqPos =  LsCum( NSeqCV ) + l
            Pknt = 0.0
            Ptot = 0.0
            Do 620 i = 1, 25, 1
               pfake(i) = 0.0
  620          continue
            p(24) = 0.0
            p(25) = 0.0
C
C ***   Adjust the raw counts for sequence residues in the analysis
C ***   group being cross validated, at the current position.
C ***   Also count the different kinds of sequence residue present.
C
            Do 650 n = 1, MaxLet, 1
               p( n ) = RawKnt( n, l )
               If( p( n ) .ge. 0.0001 )    Pknt = Pknt + MPC
  650          continue
            AAid = Seq( SeqPos )
            p( AAid ) = p( AAid ) - 1.0
            If( p( AAid ) .lt. 0.0001 )    Pknt = Pknt - MPC
C
C ***   Count the different kinds of sequence residues present
C
         Do 670 n = 1, MaxLet, 1
            If( p(n) .gt. 0.00001 )    Ptot = Ptot + 1.0
  670       continue
C
C ***   If completely ambiguous sequence character is present distribute
C ***   the appropriate counts over the other sequence characters.
C
        If( p(Unknown) .gt. 0.0001 )    Then
           Do 680 n = 1, ASize, 1
              p(n) = p(n) + p(Unknown) * BackFr(n)
  680         continue
           p(Unknown) = 0.0
        Else
        EndIf
C
C ***   If the sequences are protein adjust the counts for the ambiguous
C ***   amino acids Asx = (Asn,Asp) and Glx = (Gln,Glu).
C
         If( Protin )    Then
            p(3) = p(3) + ( p(21) / 2.0 )
            p(4) = p(4) + ( p(21) / 2.0 )
            p(21) = 0.0
            p(6) = p(6) + ( p(22) / 2.0 )
            p(7) = p(7) + ( p(22) / 2.0 )
            p(22) = 0.0
         Else
         EndIf
C
C ***   Compute the Henikoff style pseudocounts to be added to each
C ***   kind of sequence residue at this position.
C
         Do 710 n = 1, ASize, 1
            Do 700 m = 1, ASize, 1
               pfake(n) = pfake(n) + ((p(m)/Ptot) * (Qij(m,n)/QSum(m)))
  700          continue
            pfake(n) = pfake(n) * Pknt
  710       continue
C
C ***   Compute the corrected counts as the weighted average of the
C ***   actual, observed counts and the computed pseudocounts
C
         Do 730 n = 1, ASize, 1
            p(n) = (( Ptot / ( Ptot + Pknt )) * ( p(n) / Ptot ))
     *           + (( Pknt / ( Ptot + Pknt )) * ( pfake(n) / Pknt ))
  730       continue
C
C ***   Compute the entropy for the two analysis groups.  Then compute
C ***   the profile or Position Specific Scoring Matrices (log-odds matrix)
C ***   in bits.
C
            Do 750 n = 1, ASize, 1
               odds = p(n) / BackF(n,l)
C               Write( DeBugF, 10 )  n, l, odds, p(n), BackF(n,l)         DeBug
C   10 Format(/' amino acid n =', I2, ',  ungapped position l =',I4,      DeBug
C     *        ',    odds = ', F11.6,                                     DeBug
C     1       /'     corrected count p(n) =',F11.6,',    BackF(n,l) =',   DeBug
C     2         F11.6 )                                                   DeBug
               Profl1( n ) = base2 * dlog( odds )
               CVPEnt(SeqPos,1) = CVPEnt(SeqPos,1) + (p(n)*Profl1(n))
C
               odds = BackF(n,l) / p(n)
               Profl2( n ) = base2 * dlog( odds )
               CVPEnt(SeqPos,2) = CVPEnt(SeqPos,2)
     *                          + (BackF(n,l)*Profl2(n))
  750          continue
C
C ***   Compute B=Asx (21), Z=Glx (22), and X (MaxLet) averages and scores.
C
            If( Protin )    Then
               Profl1(21) = (( Profl1(3) * p(3) )
     *                    +  ( Profl1(4) * p(4) )) / ( p(3) + p(4) )
               Profl2(21) = (( Profl2(3) * p(3) )
     *                    +  ( Profl2(4) * p(4) )) / ( p(3) + p(4) )
               Profl1(22) = (( Profl1(6) * p(6) )
     *                    +  ( Profl1(7) * p(7) )) / ( p(6) + p(7) )
               Profl2(22) = (( Profl2(6) * p(6) )
     *                    +  ( Profl2(7) * p(7) )) / ( p(6) + p(7) )
            Else
            EndIf
            Profl1(MaxLet) = 0.0
            Profl2(MaxLet) = 0.0
            Max1 = 0.0
            Max2 = 0.0
            Do 800 n = 1, ASize, 1
               Profl1(MaxLet) = Profl1(MaxLet) + ( Profl1(n) * p(n) )
               Max1 = Max1 + p(n)
               Profl2(MaxLet) = Profl2(MaxLet) +(Profl2(n)*BackF(n,l))
               Max2 = Max2 + BackF(n,l)
  800          continue
            Profl1(MaxLet) = Profl1(MaxLet) / Max1
            CVPScr(SeqPos,1) = Profl1(AAid)
            Profl2(MaxLet) = Profl2(MaxLet) / Max2
            CVPScr(SeqPos,2) = Profl2(AAid)
C
            CVEntr( NSeqCV, 1 ) = CVEntr(NSeqCV,1) + CVPEnt(SeqPos,1)
            CVEntr( NSeqCV, 2 ) = CVEntr(NSeqCV,2) + CVPEnt(SeqPos,2)
            CVScor( NSeqCV, 1 ) = CVScor(NSeqCV,1) + CVPScr(SeqPos,1)
            CVScor( NSeqCV, 2 ) = CVScor(NSeqCV,2) + CVPScr(SeqPos,2)
C
 1000       continue
         If( CVScor(NSeqCV,1) .gt. thresh(DoGrup,1) )
     *       thresh(DoGrup,1) = CVScor(NSeqCV,1)
         If( CVScor(NSeqCV,2) .lt. thresh(DoGrup,2) )
     *       thresh(DoGrup,2) = CVScor(NSeqCV,2)
 2000    continue
C
      Return
      End
C
C
      Subroutine TwoEnt( NSP, TotSqP, MaxSqP, n20, NS, MxSymP, Seq,
     *                   LsCum, Qij, QSum, Pair, WithIn, Between,
     1                   AlnLen, NoGaps, LNoGap, ASize, WKnt, BKnt,
     2                   Protin, MPC, GrupID, NotID, Gap, MaxLet,
     3                   BackFr, UnKnown, DeBugF, NSPSqP )
C
      Implicit None
C
C ***   Passed Variables - Subroutine TwoEnt - calculation of the
C ***                                 entropy distance between a
C ***                                 pair of sequences.
C
      Integer       NSP, NSPSqP, TotSqP, MaxSqP, n20, NS, MxSymP,
     *              DeBugF, AlnLen, ASize, MaxLet, Gap, NotID, WKnt,
     1              BKnt, LNoGap, UnKnown
      Integer       LsCum( NSP ), Seq( TotSqP ), GrupId( NSP ),
     *              NoGaps( MaxSqP )
      Real*4        QSum( n20 ), MPC, Qij( n20, n20 ), Within( NSPSqP),
     *              Pair( NSP, NSP), Between( NSPSqP), BackFr( MxSymP )
      Logical       Protin
C
C ***   Local variables
C
      Integer       Seq1, Seq2, i, la, l, n, ni, m, aatyp1, aatyp2
      Real*8        p(25), q(25), pfake(25), qfake(25)
      Real*8        base2
C
C
      WKnt = 0
      BKnt = 0
      base2 = 1.0D 00 / dlog( 2.0D 00 )
C
      Pair( NS, NS ) = 0.0
      Do 1000 Seq1 = 1, NS-1, 1
         Pair( Seq1, Seq1 ) = 0.0
         If( GrupId( Seq1 ) .eq. 0 )    Then
            Do 110 Seq2 = Seq1+1, NS, 1
               Pair( Seq1, Seq2 ) = 0.0
               Pair( Seq2, Seq1 ) = 0.0
  110          continue
            GoTo 1000
         Else
         EndIf
         Do 900 Seq2 = Seq1+1, NS, 1 
            Pair( Seq1, Seq2 ) = 0.0
            Pair( Seq2, Seq1 ) = 0.0
            If( GrupId( Seq2 ) .eq. 0 )    GoTo 900
            Do 800 la = 1, LNoGap, 1
               l = NoGaps(la)
               aatyp1 = Seq(LsCum(seq1)+l)
               aatyp2 = Seq(LsCum(seq2)+l)
               Do 150 i = 1, Gap, 1
                  p(i) = 0.0
                  q(i) = 0.0
                  pfake(i) = 0.0
                  qfake(i) = 0.0
  150             continue
C
C ***   Check for "X" (in proteins) or "N" (in nucleic acids).  If
C ***   present set the frequencies to the average composition of the
C ***   alignment.
C
               If( aatyp1 .eq. UnKnown )    Then
                  Do 160 n = 1, ASize, 1
                     p(n) = BackFr(n)
  160                continue
               Else
                  p( aatyp1 ) = 1.0
               EndIf
               If( aatyp2 .eq. UnKnown )    Then
                  Do 170 n = 1, ASize, 1
                     q(n) = BackFr(n)
  170                continue
               Else
                  q( aatyp2 ) = 1.0
               EndIf
C
C ***   If the sequences are protein adjust the counts for the ambiguous
C ***   amino acids Asx = (Asn,Asp) and Glx = (Gln,Glu).
C
         If( Protin )    Then
            p(3) = p(3) + ( p(21) / 2.0 )
            p(4) = p(4) + ( p(21) / 2.0 )
            p(21) = 0.0
            p(6) = p(6) + ( p(22) / 2.0 )
            p(7) = p(7) + ( p(22) / 2.0 )
            p(22) = 0.0
            q(3) = q(3) + ( q(21) / 2.0 )
            q(4) = q(4) + ( q(21) / 2.0 )
            q(21) = 0.0
            q(6) = q(6) + ( q(22) / 2.0 )
            q(7) = q(7) + ( q(22) / 2.0 )
            q(22) = 0.0
         Else
         EndIf
C   10 Format('  Raw values for p(n), n = 1 -> 20;  L = ', I3 )
C   30 Format('  Raw values for p(n), n = 1 -> 20;  L = ', I3 )
C   11 Format( 5X,  10F7.1 )
C      Write( DeBugF, 10 )  l
C      Write( DeBugF, 11 ) ( p(n), n =  1, 10, 1 )
C      Write( DeBugF, 11 ) ( p(n), n = 11, 20, 1 )
C      Write( DeBugF, 30 )  l
C      Write( DeBugF, 11 ) ( q(n), n =  1, 10, 1 )
C      Write( DeBugF, 11 ) ( q(n), n = 11, 20, 1 )
C
C ***   Compute the Henikoff style pseudocounts to be added to each
C ***   kind of sequence residue at this position.
C
         Do 310 n = 1, ASize, 1
            Do 300 m = 1, ASize, 1
               pfake(n) = pfake(n) + (p(m) * (Qij(m,n)/QSum(m)))
               qfake(n) = qfake(n) + (q(m) * (Qij(m,n)/QSum(m)))
  300          continue
            pfake(n) = pfake(n) * MPC
            qfake(n) = qfake(n) * MPC
  310       continue
C    8 Format(/'  Pseudocounts for this position' )
C   12 Format('  Raw pfake(n), n = 1 -> 20;  L = ', I3 )
C   28 Format(/'  Pseudocounts for this position' )
C   32 Format('  Raw qfake(n), n = 1 -> 20;  L = ', I3 )
C   13 Format( 5X,  10F7.4 )
C      Write( DeBugF, 8 )
C      Write( DeBugF, 12 )  l
C      Write( DeBugF, 13 ) ( pfake(n), n =  1, 10, 1 )
C      Write( DeBugF, 13 ) ( pfake(n), n = 11, 20, 1 )
C      Write( DeBugF, 28 )
C      Write( DeBugF, 32 )  l
C      Write( DeBugF, 13 ) ( qfake(n), n =  1, 10, 1 )
C      Write( DeBugF, 13 ) ( qfake(n), n = 11, 20, 1 )
C
C ***   Compute the corrected counts as the weighted average of the
C ***   actual, observed counts and the computed pseudocounts
C
         Do 330 n = 1, ASize, 1
            p(n) = (( 1.0 / ( 1.0 + MPC )) * p(n) )
     *           + (( MPC / ( 1.0 + MPC )) * ( pfake(n) / MPC ))
            q(n) = (( 1.0 / ( 1.0 + MPC )) * q(n) )
     *           + (( MPC / ( 1.0 + MPC )) * ( qfake(n) / MPC ))
  330       continue
C   14 Format(/'  Divide by q= 0;   Supplemented p(n), n = 1 -> 20;',
C     *        '  L = ', I3 )
C      Write( DeBugF, 14 )  l
C      Write( DeBugF, 13 ) ( p(n), n =  1, 10, 1 )
C      Write( DeBugF, 13 ) ( p(n), n = 11, 20, 1 )
C   34 Format(/'  Divide by q= 0;   Supplemented q(n), n = 1 -> 20;',
C     *        '  L = ', I3 )
C      Write( DeBugF, 34 )  l
C      Write( DeBugF, 13 ) ( q(n), n =  1, 10, 1 )
C      Write( DeBugF, 13 ) ( q(n), n = 11, 20, 1 )
C
C ***   Compute the entropy for each of
C ***   the two analysis groups.  Then compute the two profiles or
C ***   Position Specific Scoring Matrices log-odds matrices in bits.
C ***   Also find the sequences that have the best scores in for each of
C ***   the two profiles.
C
               Do 350 n = 1, ASize, 1
C                  If( q(n) .lt. 0.0001 )    Then
C                     Write( DeBugF, 34 )  l
C                     Write( DeBugF, 13 ) ( q(ni), ni =  1, 10, 1 )
C                     Write( DeBugF, 13 ) ( q(ni), ni = 11, 20, 1 )
C                     Write( DebugF, * )  ' Seq1 =', Seq1,',   Seq2=',
C     *                                   Seq2, ',   AAtype1=', aatyp1,
C     1                                   ',   AAType2=', aatyp2,
C     2                                   ',   la=',la, ',   l=',l
C                     Stop ' Divide by q(n) = 0'
C                  Else
C                  EndIf
C                  If( p(n) .lt. 0.0001 )    Then
C                     Write( DeBugF, 14 )  l
C                     Write( DeBugF, 13 ) ( p(ni), ni =  1, 10, 1 )
C                     Write( DeBugF, 13 ) ( p(ni), ni = 11, 20, 1 )
C                     Write( DebugF, * )  ' Seq1 =', Seq1,',   Seq2=',
C     *                                   Seq2, ',   AAtype1=', aatyp1,
C     1                                   ',   AAType2=', aatyp2,
C     2                                   ',   la=',la, ',   l=',l
C                     Stop ' Divide by p(n) = 0'
C                  Else
C                  EndIf
                  Pair(Seq1,Seq2) = Pair(Seq1,Seq2)
     *                            + p(n) * base2 * dlog( p(n) / q(n) )
     1                            + q(n) * base2 * dlog( q(n) / p(n) )
  350             continue
C
C
  800          continue
            Pair(Seq2,Seq1) = Pair(Seq1,Seq2)
            If( GrupId(Seq1) .ne. NotID   .and.
     *          GrupId(Seq2) .ne. NotID )       Then
               If( GrupId(Seq1) .eq. GrupId(Seq2) )    Then
                  WKnt = WKnt + 1
                  WithIn( WKnt ) = Pair(Seq1,Seq2)
               Else
                  BKnt = BKnt + 1
                  Between( BKnt ) = Pair(Seq1,Seq2)
               EndIf
            Else
            EndIF
  900       continue
 1000    continue
C
      Return
      End
C
C
      Subroutine GetGrp( Termnl, KeyBrd, GrupFl, NS, NSP, SqName,
     *                   MxGrpP, GrupID, InvGrp, GrpNam, LGpNam,
     1                   NGroup, Title, PFile, PFLen, NotID, GrpKnt,
     2                   OutFl, DeBugF, NSqUsed, SeqUsed )
C
      Implicit None
C
C ***   Passed Variables - Subroutine GetGrp (Get Groups)
C
      Integer       Termnl, KeyBrd, GrupFl, NS, NSP, MxGrpP,
     *              NotID, GrpKnt, PFLen, OutFl, DeBugF, NSqUsed
      Integer       GrupId( NSP ), InvGrp( NSP, MxGrpP ),
     *              LGpNam( MxGrpP ), NGroup( MxGrpP ), SeqUsed( NSP )
      Character*15  SqName( NSP )
      Character*24  PFile
      Character*30  GrpNam( MxGrpP )
      Character*70  Title
      Character*80  GpFile
C
C ***   Local Variables
C
      Integer       MaxLin, LabelP, NPChrs
      Parameter  (  MaxLin = 80, NPChrs = 18, LabelP = 20 )
C
      Integer       n, k, Count, knt, ignore
      Character*1   Punct( 0:NPChrs )
      Character*2   PNumbr
      Character*8   FakeIt
      Character*15  Label( LabelP ), TLabel
      Character*30  NotUse
      Character*80  Line, UCLine
      Logical       empty, dslash
C
C ***   Function declarations
C
      Integer       LastCh
C
      Data  NotUse  /  'Sequences Ignored by Program  '  /
C
C
C ***   PUNCT contains the legal punctuation characters for VAX FORTRAN.
C ***   The first element, 0, is a sentinal needed by the binary search (POSN)
C ***   subroutine.  The punctuation characters are in ASCII collating
C ***   sequence order.  The array should be dimensioned from zero to the
C ***   number of characters which you want to allow as delimeters.
C ***   PUNCT( 0 ) = CHAR( 0 ) should always be included as part of the
C ***   program initialization as should PUNCT( 1 ) = CHAR( 9 ) for the
C ***   set of punctuation characters used here.
C
      DATA  Punct / ' ', ' ',  ' ', '!', '"', '%', '&', '''', '(',
     +             ')', '*', '+', ',', '.', '/', ':', '<',  '=',
     +             '>'                                                /
C
C
C
    1 Format(//' Output file names will all begin with the string: ',
     *          A24 )
    2 Format( A80 )
    3 Format(//' *****     Error in the file assigning sequences to',
     *         ' Groups, please correct.',
     1       //' The name of the bad file (and any prepended path) is',
     2         ' on the following line.', /' ', A80 )
    4 Format( I2 )
    5 Format( I8 )
    6 Format(//' NOTE:  You did not specify any sequences for',
     *         ' classification into',
     1        /' already defined groups.  The defined groups will be',
     2         ' cross-validated',
     3        /' and checked for consistency.', / )
    7 Format(//' Enter the name of the file describing the groups of',
     *         ' sequences from the',
     1        /' alignment just read that are to be analyzed and',
     2         ' classified.', / )
    8 Format( A24 )
    9 Format(/' Title: ', A70 )
   10 Format(/' Group Number', I4, ',     Name: ', A30,
     *        ',   Members =', I4 )
   11 Format(' Group Number',I3,',   Name Number',I4',   Name:',A15 )
   12 Format(' Seq. Number',I4,',   Group Member Number',I3,
     *       ',   Seq Name:',A15 )
   13 Format('     Sequence: ', A15 )
C
C
      Punct( 0 ) = Char( 0 )
      Punct( 1 ) = Char( 9 )
      GrpKnt = 0
      NotId = 0
      NSqUsed = 0
C
      Write( Termnl, 7 )
      Read( KeyBrd, 2 )  GpFile
      Open( Unit = GrupFl, File = GpFile, Status = 'old' )

      Do 120 n = 1, NSP, 1
         GrupId(n) = 0
         Do 110 k = 1, MxGrpP, 1
            InvGrp( n, k ) = 0
  110       continue
  120    continue
      Do 130 n = 1, MxGrpP, 1
         LGpNam( n ) = 0
         NGroup( n ) = 0
         GrpNam( n ) = ' '
  130    continue
C
C ***   Find the first non-blank line in the file and read the "problem"
C ***   name which will be used to uniquely name the output file.  The
C ***   line should begin with the keyword "problem"
C
  200 Read( GrupFl, 2 )  Line
      Call PakNam( Line, MaxLin, empty )
      If( empty )    Then
         Line(1:24) = 'GroupClassification     '
         GoTo 300
      Else
      EndIf
      UCLine = Line
      Call AllCap( UCLine, MaxLin )
  250 If( UCLine(1:7) .eq. 'PROBLEM' )    Then
         If( LGE( UCLine(8:8), 'A') .and. LLE( UCLine(8:8), 'Z' )) Then
            Line(1:7) = '       '
         Else
            Line(1:8) = '        '
         EndIf
         Call PakNam( Line, MaxLin, empty )
         If( empty )    Then
            Line(1:24) = 'GroupClassification     '
         Else
         EndIf
         PFile = Line(1:24)
         PFLen = Lastch( PFile, 24 )
         Write( Termnl, 1 )  PFile
         Write( OutFl, 1 )  PFile
      ElseIf( UCLine(1:5) .eq. 'TITLE' )    Then
         PFile = 'GroupClassification     '
         PFLen = 19
         GoTo 310
      Else
         PFile = 'GroupClassification     '
         PFLen = 19
         Title = 'Group Entropy Analysis.'
         GoTo 500
      EndIf
C
C ***   Find the title line beginnings with the keyword "title", it
C ***   must be on the first non-blank line after the "Problem" line.
C
  300 Read( GrupFl, 2 )  Line
      Call PakNam( Line, MaxLin, empty )
      If( empty )    GoTo 300
      UCLine = Line
      Call AllCap( UCLine, MaxLin )
  310 If( UCLine(1:5) .eq. 'TITLE' )    Then
         If( LGE( UCLine(6:6), 'A') .and. LLE( UCLine(6:6), 'Z' )) Then
            Line(1:5) = '     '
         Else
            Line(1:6) = '      '
         EndIf
         Call PakNam( Line, MaxLin, empty )
         Title = Line(1:70)
         Write( OutFl, 9 )  Title
      Else
         PFile = 'GroupClassification     '
         PFLen = 19
         Title = 'Group Entropy Analysis.'
         GoTo 500
      EndIf
C
C ***   Now read the sequence names for all of the groups in the data
C
  400 Read( GrupFl, 2, End = 1500 )  Line
      Call PakNam( Line, MaxLin, empty )
      If( empty )    GoTo 400
  420 UCLine = Line
      Call AllCap( UCLine, MaxLin )
  500 If( UCLine(1:5) .eq. 'GROUP' )    Then
         GrpKnt = GrpKnt + 1
         knt = 0
         If( LGE( UCLine(7:7), 'A') .and. LLE( UCLine(7:7), 'Z' )) Then
            Line(1:6) = '      '
         Else
            Line(1:7) = '       '
         EndIf
         Call PakNam( Line, MaxLin, empty )
         If( empty )    Then
            FakeIt = '        '
            Write( FakeIt, 5 )  GrpKnt
            If( GrpKnt .lt. 10 )     FakeIt(7:7) = '0'
            If( GrpKnt .lt. 100 )    FakeIt(6:6) = '0'
            FakeIt(1:5) = 'Group'
            Line(1:8) = FakeIt
         Else
         EndIf
         GrpNam(GrpKnt) = Line(1:30)
         LGpNam(GrpKnt) = LastCh( GrpNam(GrpKnt), 30 )
         UCLine = Line
         Call AllCap( UCLine, MaxLin )
         If( UCLine(1:14) .eq. 'TOBECLASSIFIED' )    NotID = GrpKnt
      Else
         GoTo 1500
      EndIf
C
C ***   get the names of sequences from the file and assign them to a group
C
  700 Read( GrupFl, 2, End = 1500 )  Line
      Call PakNam( Line, MaxLin, empty )
      If( empty )    GoTo 700
      If( Index( Line, '//' )  .gt.  0 )    Then
         dslash = .True.
      Else
         dslash = .False.
      EndIf
      Call GetLst( Termnl, KeyBrd, Count, LabelP, NPChrs, MaxLin, Line,
     *             Label, Punct )
      Do 900 k = 1, Count, 1
         TLabel = Label(k)
C         Write( DeBugF, 11 )  GrpKnt, k, TLabel
         Do 800 n = 1, NS, 1
            If( TLabel  .eq.  SqName(n) )    Then
               GrupId(n) = GrpKnt
               NGroup(GrpKnt) = NGroup(GrpKnt) + 1
               InvGrp( NGroup(GrpKnt), GrpKnt ) = n
C               Write( DeBugF, 12 )  n, NGroup(GrpKnt), SqName(n)
               GoTo 900
            Else
            EndIf
  800       continue
  900    continue
      If( .not. dslash )    Then
         GoTo 700
      ElseIf( dslash )    Then
C         Write( OutFl, 10 )  GrpKnt, GrpNam(GrpKnt), NGroup(GrpKnt)
  950    Read( GrupFl, 2, End = 1000 )  Line
         Call PakNam( Line, MaxLin, empty )
         If( empty )    GoTo 950
         GoTo 420
      EndIf
C
  
 1000 If( NotId .eq. 0 )    Write( Termnl, 6 )
C
      ignore = 0
      Do 1200 n = 1, NS, 1
        If( GrupId(n) .eq. 0 )    ignore = ignore + 1
 1200   continue
      If( ignore .ge. 1 )    Then
         Write( OutFl, 10 )  0, NotUse, ignore
         Do 1250 n = 1, NS, 1
           If( GrupId(n) .eq. 0 )    Write( OutFl, 13 )  SqName(n)
 1250      continue
      Else
      EndIf
C
      Do 1400 n = 1, GrpKnt, 1
         Write( OutFl, 10 )  n, GrpNam(n), NGroup(n)
         Do 1300 k = 1, NGroup(n), 1
            Write( OutFl, 13 )  SqName( InvGrp( k, n ) )
 1300       continue
 1400    continue
C
      Close( Unit = GrupFl, status = 'keep' )
C
C ***   Make a list of the sequences actualy used in the analysis
C
      Do 1450 n = 1, NS, 1
         If( GrupId(n) .ge. 1 )    Then
            NSqUsed = NSqUsed + 1
            SeqUsed( NSqUsed ) = n
         Else
         EndIf
 1450    continue
C
      Return
C
C ***   The following code is executed only if an end-of-file is
C ***   encountered where it should not be.
C
 1500 Write( Termnl, 3 )  GpFile
      Stop ' Bad Group Assignment File.'
      End
C
C
      SUBROUTINE GetSeq( Termnl, KeyBrd, NFile, OutFl, NSP, TotSqP,
     *                   LsTot, NS, SEQ, LSN, LsCum, Name, SeqWt,
     1                   MxSymP, Comp, Gap, MaxLet, SqPrnt, n25,
     2                   ModAA, Protin, DNA, DeBugF, AlnLen, ASize,
     3                   BackFr, UnKnown )
C
      Implicit None
C
      Integer       Termnl, KeyBrd, NFile, NSP, TotSqP, LsTot, NS,
     *              NSqTyp, OutFl, Gap, MaxLet, n25, ModAA, DeBugF,
     1              AlnLen, ASize, MxSymP, UnKnown
      Integer       Seq( TotSqP ), LSN( NSP ), Comp( 0:MxSymP ),
     *              LsCum( NSP )
      Real*4        SeqWt( NSP ), BackFr( MxSymP )
      Character*1   SqPrnt(n25), SqTyp1(2)
      Character*15  Name( NSP )
      Character*12  SeqTyp(2)
      Character*80  SqFl
      Logical       Protin, DNA
C
C ***   Local Variable Declarations
C
C ***   Note:  Raw must have the same dimensions as array Comp( 0:MxSymP )
C
      Integer       FlFmt, i, l, lb, le, lnext, n, lstart, nb, NBlock
      Integer       Raw( 0:28 ), Lets( 0:127 )
      Character*8   FlType(0:7)
      Logical       NoAlnd
C
C ***   File Type data, MUST match those in Subroutine WhatFM
C
      INTEGER      GENBNK, EMBL, NBRF, FASTA, MOLGEN, Simple, NONE,
     *             GCGMSF
C
      PARAMETER ( NONE = 0, GENBNK = 1, EMBL = 2, NBRF = 3, MOLGEN = 4,
     *            FASTA = 5, Simple = 6, GCGMSF = 7 )
C
      Data  SeqTyp  /  '  Protein   ', ' Nucleotide '  /
      Data  SqTyp1  /  'P', 'N'  /
      Data  FlType  / ' Unknown', ' GenBank', '    EMBL', 'PIR/NBRF',
     *                '  Molgen', '   Fasta', '  Simple', ' GCG MSF' /
C
C ***   The Lets array translates the sequences from a letter to a numeric
C ***   representation.  Each letter is replaced by its location in the
C ***   alphabetic collating sequence.
C
C                            ' '       #           *         -   .
      DATA  Lets  / 32 * -1,  0, 0, 0, 0,3*-1,3*0,100, 0,0, 27, 27, 0,     HBN
C
C                        1   2                           A   B   C   D
     *             -1, 201,102, 10 * -1, 0, -1, -1, -1,  1,  2,  3,  4,    HBN
C
C                   E    F   G   H   I   J   K   L   M   N   O   P   Q
     1              5,   6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,    HBN
C
C                   R    S    T   U   V   W   X   Y   Z             a
     2             18,  19,  20, 21, 22, 23, 24, 25, 26, 6 *  0,    1,     HBN
C
C                   b   c   d   e   f   g   h   i   j   k   l   m   n
     3              2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,     HBN
C
C                   o   p   q   r   s   t   u   v   w   x   y   z
     4             15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,         HBN
     5              5 *  0   /                                             HBN
C
C
    1 Format(//' Enter the name of the file holding the sequences.'/)
    2 Format(/' The data consist of', I5, A12, 'sequences.'/)
    3 Format( A80 )
    4 Format(/' Your sequence data file:', /' ', A80,
     *       /' is in', A8, ' format.')
    5 Format( ' ', A10, '     Length =', I6 )
    6 Format( '    ', 25I3: )
    7 Format( '    ', 7( 1X,10A1: ) )
    8 Format(/'  TestAlign.msf    MSF: ', I5,'   Type: ',A1,
     *        '    Check: ', '  3679    ..', // )
    9 Format('  Name: ', A12, '   Len: ', I4, '   Check: ', I4,
     *       '   Weight: ', F7.4 )
   10 Format(/, '//' )
   11 Format(/' ', 10X, I6, 50X, I4 )
   12 Format( A12, 1X, 5(1X, 10A1:) )
   13 Format(' ')
C
C
      Write( Termnl, 1 )
      Read( KeyBrd, 3 )  SqFl
      Open( Unit = NFile, File = SqFl, Status = 'old' )
C
      Call WhatFm( NFile, Termnl, FlFmt )
      Write( OutFl, 4 )  SqFl, FlType(FlFmt)
C
      If( FlFmt .eq. GCGMSF )    Then
            Call RdMSF( Termnl, NFile, NSP, TotSqP, LsTot, NS, Seq,
     *                  LSN, Name, MxSymP, Raw, Lets, KeyBrd, AlnLen,
     1                  LsCum, SeqWt )
      ElseIf( FlFmt .eq. GenBnk )    Then
            Call ReadGB( TERMNL, NFILE, NSP, TotSqP, LsTot, NS, SEQ,
     *                   LSN, Name, MxSymP, Raw, Lets)
      ElseIf( FlFmt .eq. EMBL )    Then
            Call RdEMBL( TERMNL, NFILE, NSP, TotSqP, LsTot, NS, SEQ,
     *                   LSN, Name, MxSymP, Raw, Lets)
      ElseIf( FlFmt .eq. NBRF )    Then
            Call RdPIR( TERMNL, NFILE, NSP, TotSqP, LsTot, NS, SEQ,
     *                  LSN, Name, MxSymP, Raw, Lets)
      ElseIf( FlFmt .eq. FASTA )    Then
            Call RFasta( TERMNL, NFILE, NSP, TotSqP, LsTot, NS, SEQ,
     *                   LSN, Name, MxSymP, Raw, Lets)
      ElseIf( FlFmt .eq. MolGen )    Then
            Call MolGn( TERMNL, NFILE, NSP, TotSqP, LsTot, NS, SEQ,
     *                  LSN, Name, MxSymP, Raw, Lets)
      ElseIf( FlFmt .eq. Simple )    Then
            Call Simpl( TERMNL, NFILE, NSP, TotSqP, LsTot, NS, SEQ,
     *                  LSN, Name, MxSymP, Raw, Lets)
      ElseIf( FlFmt .eq. None )    Then
      EndIf
      Close( Unit = NFile, Status = 'keep' )
C
      Call FixSeq( LsTot, TotSqP, Seq, MxSymP, Raw, Comp, Gap, MaxLet,
     *             ASize, SqPrnt, n25, ModAA, Protin, DNA, BackFr,
     1             UnKnown )
C
      If( Protin )    Then
         NSqTyp = 1
      Else
         NSqTyp = 2
      EndIf
      Write( OutFl, 2 )  NS, SeqTyp(NSqTyp)
      NoAlnd = .False.
      If( FlFmt .ne. GCGMSF )    Then
         AlnLen = LSN(1)
         LsCum(1) = 0
         SeqWt(1) = 1.0
         Do 100 n = 2, NS
            LsCum(n) = LsCum(n-1) + Lsn(n-1)
            SeqWt(n) = 1.0
            If( LSN(n) .eq. AlnLen )    Then
            ElseIf( LSN(n) .gt. AlnLen )    Then
               NoAlnd = .True.
               AlnLen = LSN(n)
            ElseIf( LSN(n) .lt. AlnLen )    Then
               NoAlnd = .True.
            EndIf
  100       continue
      Else
      EndIf
C
C ***   Debugging code to check the correct translation of sequences
C
C      lnext = 0
C      Do 210 n = 1, NS
C         Write( DeBugF, 5 )  Name(n), Lsn(n)
C         lb = -24
C  200    lb = lb + 25
C         le = lb + 24
C         If( le .gt. Lsn(n) )    le = Lsn(n)
C         Write( DeBugF, 6 )  ( Seq(lnext+l), l = lb, le, 1 )
C         If( le .lt. Lsn(n) )    GoTo 200
C         lnext = lnext + Lsn(n)
C  210    continue
C
C***   Write back translation of sequences to results file
C
      If( NoAlnd )    Then
         lnext = 0
         Do 250 n = 1, NS
            Write( OutFl, 5 )  Name(n), Lsn(n)
            lb = -69
  230       lb = lb + 70
            le = lb + 69
            If( le .gt. Lsn(n) )    le = Lsn(n)
            Write( OutFl, 7 )  ( SqPrnt(Seq(lnext+l)), l = lb, le, 1 )
            If( le .lt. Lsn(n) )    GoTo 230
            lnext = lnext + Lsn(n)
  250       continue
      Else
         Write( OutFl, 8 )  AlnLen, SqTyp1(NSqTyp)
         Do 310 n = 1, NS, 1
            Write( OutFl, 9 ) Name(n), LSN(n), n, SeqWt(n)
  310       continue
         Write( OutFl, 10 )
         NBlock = ( AlnLen + 49 ) / 50
         lb = -49
         Do 350 nb = 1, NBlock, 1
            lb = lb + 50
            le = lb + 49
            If( le .gt. AlnLen )    le = AlnLen
            Write( OutFl, 11 )  lb, le
            Do 340 n = 1, NS, 1
               i = LsCum(n)
               Write(OutFl,12)  Name(n), ( SqPrnt(Seq(i+l)),l=lb,le,1)
  340          continue
  350       continue
         Write( OutFl, 13 )
      EndIf
C
      Return
      End



C
      Subroutine RdMSF( Termnl, MSFfl, NSeqP, SqBufP, LSqTot, NSeq,
     *                  Seq, SqLen, SqName, MxSymP, Raw, Lets, KeyBrd,
     1                  AlnLen, OffSet, SeqWt )
C
      Implicit None
C
C ***   Passed Variable and Parameter declarations
C
      Integer       Termnl, KeyBrd, MSFfl, NSeqP, SqBufP, LSqTot, NSeq,
     *              AlnLen, MxSymP
      Integer       Seq( SqBufP ), SqLen( NSeqP ), OffSet( NSeqP ),
     *              Raw( 0:MxSymP ), Lets( 0:127 )
      Real*4        SeqWt( NSeqP )
      Character*15  SqName( NSeqP )
C
C ***   Local Variable and Parameter Declarations
C
C
C ***   Variables used in parsing input lines
C ***   NSeqP2 should be set to at least as large as NSeqP that is
C ***   passed in from the main program
C
      INTEGER       NumbrP, MaxLin, NSeqP2
C
      PARAMETER  (  NumbrP = 20,   MaxLin = 80, NSeqP2 = 600  )
C
      INTEGER       Count
      INTEGER       Digits( NumbrP ), LayOut( NumbrP )
      CHARACTER*31  TmpStr, Numbrs( NumbrP )
      CHARACTER*80  ILine                   !   must be character*MaxLin
C
C ***   Variable declarations used in converting character strings to their
C ***   numeric values
C
      INTEGER     FmtLen
C
      PARAMETER ( FmtLen = 20 )
C
      INTEGER        DPoint, MaxLen, SeqInt, knt
      CHARACTER*20   FmtVar                  !   must be character*FmtLen
C
      INTEGER       Len, I, n, nb, ns, NBlock
      Integer       LnSeq( NSeqP ), SqIDLn( NSeqP2 )
      Logical       PrSeq, SqType, BadSeq
C
C ***   functions and subroutines
C
      INTEGER      LastCh
C
      EXTERNAL     NUMLST, ALLCAP
C

C
    2 FORMAT ( A80 )
    3 Format(//' The length specification in the following data line',
     *         ' is not an integer.',
     1        /' Please run the GCG Reformat program to correctly',
     2         ' format your .msf file.', / )
    4 Format(//' The weight specification in the following data line',
     *         ' is not a decimal',
     1        /' number.  Please run the GCG Reformat program to',
     2         ' correctly format your .msf',
     3        /' file.', / )
    5 Format(//' The weight specification in the following data line',
     *         ' is missing or',
     1        /' incorrect.  Please run the GCG Reformat program to',
     2         ' correctly format your .msf',
     3        /' file.', / )
    6 Format(/)
    7 Format(//' The name, ', A12, ', at the beginning of a line of',
     *         ' sequence data did not',
     1        /' match any of the names read from the data section of',
     2         ' the .msf file.', / )
    8 Format(//' The length field, "MSF: length", in the header line',
     *         ' is not an integer.',
     1        /' Please run the GCG Reformat program to correctly',
     2         ' format your .msf file.', / )
    9 Format(//' *****  ERROR  *****  The sequence ', A12, ' has',
     *         ' fewer characters (', I5, ')',
     1        /' than were expected (', I5, ') from the length',
     2         ' specification in the file header.', / )
   10 Format(//' *****  ERROR  *****  The sequence ', A12, ' has',
     *         ' more characters (', I5, ')',
     1        /' than were expected (', I5, ') from the length',
     2         ' specification in the file header.', / )
   11 Format(//' *****  ERROR  *****  The expected line terminated by',
     *         ' double dots ".." that',
     1        /' marks the end of the comment section and the',
     2         ' beginning of the data section',
     3        /' was not found.  Please make sure that you specifed',
     4         ' a GCG .msf formatted file.', / )
   12 Format(//' *****  ERROR  *****  The expected line terminated by',
     *         ' double slash, "//"',
     1        /' that marks the end of the header section and the',
     2         ' beginning of the sequence',
     3        /' data section was not found.  Please make sure that',
     4         ' you specifed a GCG .msf',
     5        /' formatted file.', / )
C   13 Format('  Name: ', A12, ',    Length: ', I5, ',   Weight: ',        DeBug
C     *          F7.4 )                                                    DeBug
C   14 Format('  Alignment Length: ', I8 )                                 DeBug
   15 Format(' ', A70 )
C   16 Format(' Block:',I4,'    Line:',I4,',  Data: ', A40 )               DeBug
C
C
      Do 100 i = 0, MxSymP, 1
         Raw(i) = 0
  100    continue
      LSqTot = 0
C
C ***   Find the line marking the beginning of the data.  This line ends
C ***   with "double dots, .."
C
  200 Read( MSFfl, 2, End = 2000 )  ILine
         If( Index( ILine, 'MSF: ' )   .le. 0    .or.
     *       Index( ILine, 'Check: ' ) .le. 0    .or.
     *       Index( ILine, '..' ) .le.  0  )    Then
            GoTo 200
         Else
            Call NumLst( Termnl, KeyBrd, Count, NumbrP, MaxLin, ILine,
     *                   Numbrs, Digits, LayOut )
            Do 250 n = 1, Count, 1
               If( Numbrs(n)(1:Digits(n)) .eq. 'MSF:' )    Then
                  If( LayOut( n+1 )   .eq.   1 )    Then        !   Integer
                     Call IForm( Termnl, Digits(n+1), FmtLen, FmtVar )
                     Read( Numbrs(n+1), FmtVar )  AlnLen
                  Else
                     Write( Termnl, 8 )
                     Len = LastCh( ILine, 80 )
                     Call AForm( Termnl, Len, FmtLen, FmtVar )
                     Write( Termnl, FmtVar )  ILine(1:len)
                     Stop ' Error in alignment length field.'
                  EndIf
               ElseIf( Numbrs(n)(1:Digits(n)) .eq. 'Type:' )    Then
                  If( Numbrs(n+1) .eq. 'P'   .or.
     *                Numbrs(n+1) .eq. 'p' )    Then
                     PrSeq = .True.
                  Else
                     PrSeq = .False.
                  EndIf
                  SqType = .True.
                  GoTo 300
               Else
               EndIf
  250          continue
            SqType = .False.
            PrSeq = .False.
         EndIf
C
C ***   Now get Names, Lengths, and Weights for each sequence
C
  300 NSeq = 0
C      Write( Termnl, 14 )  AlnLen                                         DeBug
  310 Read( MSFfl, 2, End = 2010 )  ILine
         If( Index( ILine, 'Name:' ) .gt. 0 )    Then
            Call NumLst( Termnl, KeyBrd, Count, NumbrP, MaxLin, ILine,
     *                   Numbrs, Digits, LayOut )
C
C ***   Get the sequence identification information
C
            NSeq = NSeq + 1
            Do 350 n = 1, Count, 1
               TmpStr = Numbrs(n)
               Call AllCap( TmpStr, Digits(n) )
               If( TmpStr(1:5) .eq. 'NAME:' )    Then
                  SqName(NSeq) = Numbrs(n+1)(1:12)
                  SqIDLn(NSeq) = Digits(n)
                  If( SqIDLn(NSeq) .gt. 12 )    SqIDLn(NSeq) = 12
               ElseIf( TmpStr(1:4) .eq. 'LEN:' )    Then
                  If( LayOut( n+1 )   .eq.   1 )    Then        !   Integer
                     Call IForm( Termnl, Digits(n+1), FmtLen, FmtVar )
                     Read( Numbrs(n+1), FmtVar )  SqLen(NSeq)
                  Else
                     Write( Termnl, 3 )
                     Len = LastCh( ILine, 80 )
                     Call AForm( Termnl, Len, FmtLen, FmtVar )
                     Write( Termnl, FmtVar )  ILine(1:len)
                     Stop ' Error in sequence length data.'
                  EndIf
               ElseIf( TmpStr(1:7) .eq. 'WEIGHT:' )    Then
                  If( LayOut(n+1)   .EQ.   2 )    Then    !   F formatted real
                     DPoint = Index( Numbrs( n+1 ), '.' )
                     Call FForm( Termnl, Digits(n+1), DPoint, FmtLen,
     *                           FmtVar )
                     Read( Numbrs(n+1), FmtVar )  SeqWt(NSeq)
                     GoTo 310
                  Else
                     Write( Termnl, 4 )
                     Len = LastCh( ILine, 80 )
                     Call AForm( Termnl, Len, FmtLen, FmtVar )
                     Write( Termnl, FmtVar )  ILine(1:len)
                     Stop ' Error in sequence weight data.'
                  EndIf
               Else
               EndIf
  350          continue
            Write( Termnl, 5 )
            Len = LastCh( ILine, 80 )
            Call AForm( Termnl, Len, FmtLen, FmtVar )
            Write( Termnl, FmtVar )  ILine(1:len)
            Stop ' Error - missing sequence weight.'
         ElseIf( ILine(1:2) .eq. '//' )    Then
            GoTo 400
         Else
            GoTo 310
         EndIf
C
C ***   Now read the sequence data, first find the maximum sequence
C ***   length and set all buffer regions to the length of the sequence.
C
  400 MaxLen = 0
C      Do 401 ns = 1, NSeq, 1                                              DeBug
C         Write( Termnl, 13 )  SqName(ns), SqLen(ns), SeqWt(ns)            DeBug
C  401    continue                                                         DeBug
      Do 410 ns = 1, NSeq, 1
         LnSeq(ns) = 0
         If( SqLen(ns) .gt. MaxLen )    MaxLen = SqLen(ns)
  410    continue
      AlnLen = MaxLen
      OffSet(1) = 0
      Do 420 ns = 2, NSeq, 1
         OffSet(ns) = Offset(ns-1) + SqLen(ns-1)
  420    continue
C
C ***   Now read the sequence data
C
      NBlock = ( MaxLen + 49 ) / 50
      Do 490 nb = 1, NBlock, 1
         Read( MSFfl, 6 )
         Do 480 ns = 1, NSeq, 1
            Read( MSFfl, 2 )  ILine
C            Write( Termnl, 16 )  nb, ns, ILine(1:40)                      DeBug
            Call NumLst( Termnl, KeyBrd, Count, NumbrP, MaxLin, ILine,
     *                   Numbrs, Digits, LayOut )
            If( Numbrs(1)(1:12) .ne. SqName(ns) )    Then
               Write( Termnl, 7 )  Numbrs(1)(1:12)
               Write( Termnl, 15 )  ILine(1:70)
               Stop ' Unrecognized Sequence Name.'
            Else
            EndIf
            Do 450 knt = 2, Count, 1
               Do 440 i = 1, Digits(knt), 1
                  SeqInt = Lets( Ichar( Numbrs(knt)(i:i) ) )
                  If( SeqInt .ge. 1   .and.   SeqInt .le. MxSymP ) Then
                     LnSeq(ns) = LnSeq(ns) + 1
                     Seq( OffSet(ns) + LnSeq(ns) ) = SeqInt
                     Raw( SeqInt ) = Raw( SeqInt ) + 1
                     LSqTot = LSqTot + 1
                  Else
                  EndIf
  440             continue
  450          continue
  480       continue
  490    continue
C
C ***   Now check to see if the actual sequence lengths match the lengths
C ***   writen in the data section of the file.
C
      BadSeq = .False.
      Do 550 ns = 1, NSeq, 1
         If( LnSeq(ns) .eq. SqLen(ns) )    Then
         ElseIf( LnSeq(ns) .lt. SqLen(ns) )    Then
            Write( Termnl, 9 )  SqName(ns), LnSeq(ns), SqLen(ns)
            BadSeq = .True.
         ElseIf( LnSeq(ns) .gt. SqLen(ns) )    Then
            Write( Termnl, 10 )  SqName(ns), LnSeq(ns), SqLen(ns)
            BadSeq = .True.
         EndIf
  550    continue
      If( BadSeq )    Stop ' Corrupted .msf file, please reformat.'
C
      Return
 2000 Write( Termnl, 11 )
      Stop ' File may not be GCG .msf format.'
 2010 Write( Termnl, 12 )
      Stop ' File may not be GCG .msf format.'
      END
C
C
      SUBROUTINE NUMLST( Termnl, KeyBrd, Count, NumbrP, MaxLin, ILine,
     *                   Numbrs, Digits, LayOut )
C
C ****************************************************************************
C **   Copyright  --  Pittsburgh Supercomputing Center  --  December 1988   **
C ********************  Written by Hugh B. Nicholas  *************************
C
C ***   Subroutine NUMLST gets a list of tokens from a line of input, ILine,
C ***   and stores the individual tokens in an array Numbrs.  If the last
C ***   character in the line is < another line of is solicited to include
C ***   in the list (Numbrs).  A token is defined as a string of characters
C ***   delimited by spaces or commas.  Integer function STYLE classifies
C ***   the tokens as integer or real numbers or as an alphanumeric string.
C ***   This result is returned in array LayOut.  The maximum number of
C ***   characters allowed on a line is set by MaxLin.  Count returns the
C ***   number of tokens found.  Digits hold the number of characters in a
C ***   token.  The parameter MAXCHR defines the maximum number of characters
C ***   allowed in a token.  If MAXCHR is redefined, the the character
C ***   variable NUMBER should be redeclared to be the same length.  The
C ***   passed array Numbrs should be declared to have elements of MAXCHR
C ***   characters.
C
      IMPLICIT NONE
C
      INTEGER       MAXCHR
C
      PARAMETER  (  MAXCHR = 31  )
C
      INTEGER       Termnl, KeyBrd, Count, MaxLin, NumbrP
      INTEGER       Digits( NumbrP ), LayOut( NumbrP )
      CHARACTER*(*) Numbrs( NumbrP )
      CHARACTER*(*) ILine
C
      INTEGER       I, K, IBGN, IEND, NEXTC, LASTC, LENGTH, LSCAN
C      CHARACTER*10  FMT2
      CHARACTER*31  NUMBER
      LOGICAL       EMPTY, COMMA, DEFALT
C
      INTEGER   LASTCH, STYLE                          !   function names
C
C
C      DATA  FMT2  /  '(A       )'  /            !   3 : 9 out of 10
C
C
    1 FORMAT ( //' Enter an additional line to be processed.  Enter',
     *         ' an < at the end of',
     1        /' the line if you need to continue the input onto',
     2         ' another line.',//)
C
C    3 FORMAT ( //'     **********     WARNING     **********',
C     *          /' A completely blank line was entered - please',
C     1           ' respond with valid input.'
C     2          /' A second blank line will cause the program to',
C     3           ' ABORT.'//)
C
C
C      CALL AFIELD( Termnl, MaxLin, 7, FMT2(3:9) )
   90 Count = 0
      NEXTC = 1
      COMMA = .TRUE.
      LASTC = LASTCH( ILine, MaxLin )
C
  100 CALL GETNUM( NEXTC, LASTC, ILine, IBGN, IEND, LSCAN, EMPTY, COMMA,
     *             DEFALT, NUMBER, LENGTH, MAXCHR )
C
      IF ( .NOT. EMPTY )    THEN
C
C ***       Another "token" has been returned from subroutine getnum, add
C ***       it to the list of possibly numeric "tokens" found on the 
C ***       input line.
C
         Count = Count + 1
         Numbrs( Count ) = NUMBER
         Digits( Count ) = LENGTH
         IF( .NOT. DEFALT )    THEN
            LayOut( Count ) = STYLE( NUMBER, LENGTH )
         ELSE
            LayOut( Count ) = 0
         END IF
         NEXTC = LSCAN + 1
         IF( LSCAN .LT. LASTC )    GO TO 100
C
C      ELSE IF( EMPTY   .AND.   Count .GT. 0 )    THEN
C
C         IF( ILine(LASTC:LASTC) .EQ. '<' )    THEN   !  continuation requested
C            Write( Termnl, 1 )
C            Read( KeyBrd, FMT2 )  ILine
C            NEXTC = 1
C            LASTC = LASTCH( ILine, MaxLin )
C            GO TO 100
C         ELSE
C         END IF
C
      ELSE IF( EMPTY   .AND.   Count .LE. 0 )    THEN
C
         RETURN
C
C ***   The code below is appropriate for interactive mode use - it responds
C ***   to a blank line of input by reading another line and stopping the
C ***   program if the second line is also blank.  It is an alternative to
C ***   the simple RETURN statement above.  If you use this also removed the
C ***   "C's" inactivating format 3 above.
C
C         Write( Termnl, 3 )
C         Read( KeyBrd, FMT2 )  ILine
C         LASTC = LASTCH( ILine, MaxLin )
C         IF( LASTC .EQ. 0 )   STOP ' invalid input.'
C         GO TO 90
C
      END IF
C
      RETURN
      END
C
         SUBROUTINE GETNUM( BGN, LAST, LINE, START, FINISH, LSCAN,
     *                      EMPTY, COMMA, DEFALT, NUMBER, LENGTH,
     1                      LIMIT )
C
C ****************************************************************************
C **   Copyright  --  Pittsburgh Supercomputing Center  --  December 1988   **
C ********************  Written by Hugh B. Nicholas  *************************
C
C
C ***   Subroutine GETNUM gets the next token from a character string.  The
C ***   subroutine begins searching at position BGN of the string and stops
C ***   at position LAST.  If no token is found in this region of the line
C ***   the logical variable EMPTY is set .TRUE..  If a token is found, the
C ***   variable EMPTY is set .FALSE. and the token copied to the variable
C ***   NUMBER for return to the calling program.  If a token is found, the
C ***   variables START and FINISH hold the positions of the first and last
C ***   characters of the token in the string (LINE).  Tokens longer than 
C ***   LIMIT characters are truncated to LIMIT characters.  The number
C ***   of characters in the token is returned as LENGTH.  The logical
C ***   variable COMMA is returned as .TRUE. if the delimiter following a
C ***   number was a comma, ','.  The logical variable DEFALT is returned as
C ***   .TRUE. if there are two commas without an intervening token in a
C ***   string.  LSCAN records the last position in the string that has been
C ***   examined.
C
      IMPLICIT NONE
C
      INTEGER        BGN, LAST, START, FINISH, LENGTH, LSCAN, LIMIT
      CHARACTER*(*)  NUMBER
      CHARACTER*(*)  LINE
      LOGICAL        EMPTY, COMMA, DEFALT
C
      INTEGER   I, LIMITM
      INTEGER   POSN                                       !   function name
C
C
      LIMITM = LIMIT - 1
      START = 0
      FINISH = 0
      DEFALT = .FALSE.
      NUMBER = '                               '
C
      IF( BGN .GT. LAST )    THEN
         EMPTY = .TRUE.
         RETURN
      ELSE
      END IF
C
C ***   Find the first non-punctuation character on the line
C
      DO 120 I = BGN, LAST
         IF(  LINE(I:I) .NE. ' '   .AND.
     *        LINE(I:I) .NE. ','   .AND.
     1        LINE(I:I) .NE. '<'  )      THEN
            START = I
            EMPTY = .FALSE.
            COMMA = .FALSE.
            GO TO 150
         ELSE IF( LINE(I:I) .EQ. ',' )    THEN
            IF( COMMA )    THEN
               DEFALT = .TRUE.
               EMPTY = .FALSE.
               START = I
               FINISH = I - 1
               LSCAN = I
               LENGTH = 0
               RETURN                              !   with number = spaces
            ELSE IF( .NOT. COMMA )    THEN
               COMMA = .TRUE.
            END IF
         ELSE                 !   blank space
         END IF
  120    CONTINUE
      EMPTY = .TRUE.
      RETURN
C
C ***   Find the first punctuation character after the START location.
C
  150 DO 200 I = START + 1, LAST
         IF( LINE(I:I) .EQ. ' ' )    THEN
            LSCAN = I
            FINISH = I - 1
            LENGTH = FINISH - START + 1
            NUMBER(1:LENGTH) = LINE(START:FINISH)
            RETURN
         ELSE IF( LINE(I:I) .EQ. ',' )    THEN
            COMMA = .TRUE.
            LSCAN = I
            FINISH = I - 1
            LENGTH = FINISH - START + 1
            NUMBER(1:LENGTH) = LINE(START:FINISH)
            RETURN
         ELSE
         END IF
  200    CONTINUE
C
      LSCAN = LAST
      IF( LAST - START  .LE. LIMITM )    THEN
         FINISH = LAST
      ELSE
         FINISH = START + LIMITM
      END IF
      LENGTH = FINISH - START + 1
      NUMBER(1:LENGTH) = LINE(START:FINISH)
C
      RETURN
      END
C
      SUBROUTINE AFORM( Termnl, NCHRS, MAX, FmtVar )
C
C ****************************************************************************
C **   Copyright  --  Pittsburgh Supercomputing Center  --  December 1988   **
C ********************  Written by Hugh B. Nicholas  *************************
C
C
C ***   Subroutine AFORM fills a character variable with a format
C ***   specification for A, I, F, E, or D editing descriptors.  Each
C ***   different type of edit descriptor is associated with its own
C ***   entry point.
C
C     ***************     Variable Definitions     ***************
C
C ***  DPoint = The location of the decimal point in the character string
C               representation of a floating point number.
C ***  FmtVar = The character variable to hold the format specification.
C ***  LPos   = The position of the letter E or D in the chartacter string
C               representation of an exponential number.
C ***  NCHRS  = The total number of characters to be read with the format.
C ***  MAX    = The number of spaces (bytes) in the format variable (FmtVar).
C ***  MOST   = The largest number of significant Digits that can fit into
C               a given field size (NCHRS) in an E or D format specification.
C
C
      Implicit None
C
      INTEGER         Termnl, NCHRS, MAX, DPoint, LPos
      CHARACTER*(*)   FmtVar
C
C ***   Local Variables
C
      INTEGER         I, Count, DCount, DF, MOST, POS
      CHARACTER*15    SIZE, DSIZE
      LOGICAL         EMPTY
C
      INTEGER         LASTCH                        !   function declaration
C
C ***   The following EXTERNAL statement list all of the user supplied
C ***   subroutines and entry points called in subroutine AFORM.
C
      EXTERNAL  PAKNAM
C
C
    1 FORMAT ( I15 )
    2 FORMAT ( /' ', I10, ' Characters were asked for in a format but',
     *         /' there was only space for', I10 )
C
C
      Write( SIZE, 1 )  NCHRS
      CALL PAKNAM( SIZE, 15, EMPTY )
      Count = LASTCH( SIZE, 15 )
      FmtVar(1:2) = '(A'
      IF( Count .LE. MAX - 3 )    THEN
         FmtVar(3:Count+2) = SIZE(1:Count)
         FmtVar(Count+3:Count+3) = ')'
         DO 100 I = Count + 4, MAX, 1
  100       FmtVar(I:I) = ' '
      ELSE
         Write(Termnl,2)  Count, MAX
         DO 110 I = 3, MAX - 1
  110       FmtVar(I:I) = '9'
         FmtVar(MAX:MAX) = ')'
      END IF
      RETURN
C
C
      ENTRY  IFORM( Termnl, NCHRS, MAX, FmtVar )
C
C
      Write( SIZE, 1 )  NCHRS
      CALL PAKNAM( SIZE, 15, EMPTY )
      Count = LASTCH( SIZE, 15 )
      FmtVar(1:2) = '(I'
      IF( Count .LE. MAX - 3 )    THEN
         FmtVar(3:Count+2) = SIZE(1:Count)
         FmtVar(Count+3:Count+3) = ')'
         DO 150 I = Count + 4, MAX, 1
  150       FmtVar(I:I) = ' '
      ELSE
         Write(Termnl,2)  Count, MAX
         DO 160 I = 3, MAX - 1
  160       FmtVar(I:I) = '9'
         FmtVar(MAX:MAX) = ')'
      END IF
      RETURN
C
C
      ENTRY  FFORM( Termnl, NCHRS, DPoint, MAX, FmtVar )
C
C
      Write( SIZE, 1 )  NCHRS
      CALL PAKNAM( SIZE, 15, EMPTY )
      Count = LASTCH( SIZE, 15 )
      DF = NCHRS - DPoint
      Write( DSIZE, 1 )  DF
      CALL PAKNAM( DSIZE, 15, EMPTY )
      DCount = LASTCH( DSIZE, 15 )
      FmtVar(1:2) = '(F'
      IF( Count + DCount  .LE.  MAX - 4 )    THEN
         FmtVar(3:Count+2) = SIZE(1:Count)
         POS = Count + 3
         FmtVar(POS:POS) = '.'
         FmtVar(POS+1:POS+DCount) = DSIZE(1:DCount)
         POS = POS + DCount + 1
         FmtVar(POS:POS) = ')'
         DO 200 I = POS + 1, MAX, 1
  200       FmtVar(I:I) = ' '
      ELSE
         Write(Termnl,2)  Count, MAX
         DO 210 I = 3, MAX - 3
  210       FmtVar(I:I) = '9'
         FmtVar(MAX-2:MAX) = '.0)'
      END IF
      RETURN
C
C
      ENTRY  EFORM( Termnl, NCHRS, LPos, MAX, FmtVar )
C
C
      Write( SIZE, 1 )  NCHRS
      CALL PAKNAM( SIZE, 15, EMPTY )
      Count = LASTCH( SIZE, 15 )
      DF = LPos
      MOST = NCHRS - 7
      IF( DF .GT. MOST )    DF = MOST
      Write( DSIZE, 1 )  DF
      CALL PAKNAM( DSIZE, 15, EMPTY )
      DCount = LASTCH( DSIZE, 15 )
      FmtVar(1:2) = '(E'
      IF( Count + DCount  .LE.  MAX - 4 )    THEN
         FmtVar(3:Count+2) = SIZE(1:Count)
         POS = Count + 3
         FmtVar(POS:POS) = '.'
         FmtVar(POS+1:POS+DCount) = DSIZE(1:DCount)
         POS = POS + DCount + 1
         FmtVar(POS:POS) = ')'
         DO 250 I = POS + 1, MAX, 1
  250       FmtVar(I:I) = ' '
      ELSE
         Write(Termnl,2)  Count, MAX
         DO 260 I = 3, MAX - 3
  260       FmtVar(I:I) = '9'
         FmtVar(MAX-2:MAX) = '.0)'
      END IF
      RETURN
C
C
      ENTRY  DFORM( Termnl, NCHRS, LPos, MAX, FmtVar )
C
C
      Write( SIZE, 1 )  NCHRS
      CALL PAKNAM( SIZE, 15, EMPTY )
      Count = LASTCH( SIZE, 15 )
      DF = LPos
      MOST = NCHRS - 7
      IF( DF .GT. MOST )    DF = MOST
      Write( DSIZE, 1 )  DF
      CALL PAKNAM( DSIZE, 15, EMPTY )
      DCount = LASTCH( DSIZE, 15 )
      FmtVar(1:2) = '(D'
      IF( Count + DCount  .LE.  MAX - 4 )    THEN
         FmtVar(3:Count+2) = SIZE(1:Count)
         POS = Count + 3
         FmtVar(POS:POS) = '.'
         FmtVar(POS+1:POS+DCount) = DSIZE(1:DCount)
         POS = POS + DCount + 1
         FmtVar(POS:POS) = ')'
         DO 300 I = POS + 1, MAX, 1
  300       FmtVar(I:I) = ' '
      ELSE
         Write(Termnl,2)  Count, MAX
         DO 310 I = 3, MAX - 3
  310       FmtVar(I:I) = '9'
         FmtVar(MAX-2:MAX) = '.0)'
      END IF
C
C
      RETURN
      END
C
      SUBROUTINE READGB( TERMNL, NFILE, NSP, TotSqP, LsTot, NS, SEQ,
     *                   LSN, Name, MxSymP, Raw, Lets)
C
C ***   Subroutine READGB reads a genbank formatted ASCII sequence file
C       and returns the nucleic acid sequences, coded in numeric form,
C       stored contiguously in array SEQ.  The lengths of the sequences
C       are returned in array LSN.  The GenBank locus names are returned
C       in the character array Name 
C
C       TERMNL is the FORTRAN unit number of the standard output file, where
C       error messages should be sent.  NFILE is the FORTRAN unit number for
C       the GenBank sequence data file.
C
      IMPLICIT NONE
C
C ***   Passed Variable declarations
C
      INTEGER       TERMNL, NFILE, NSP, TotSQP, LsTot, NS, MxSymP
      INTEGER       SEQ( TotSQP ), Raw ( 0:MxSymP ), LSN( NSP ),
     *              Lets( 0:127 )
      CHARACTER*15  Name( NSP )
C
C ***   Local variable declarations
C
      INTEGER       FIRST, I, LC, LS, NA, NDEND, ND, NP, jbgn, jend,
     *              j, jb, jblock, Length
      CHARACTER*1   LSEQ( 60 )
      CHARACTER*6   ACNUM
      CHARACTER*80  SEQDEF, LOCLN, ACC, ORIGIN, TERM, DEFLN
C
C
    1 FORMAT ( A80 )
    2 FORMAT ( 1 ( 9X, 6 ( 1X, 10( A1: ) ) ) )
    3 FORMAT ( //' An ACCESSION line was not found for the entry:',
     *         2 ( /' ', A80 ), / )
    4 FORMAT (//' No TERNIMATOR "//" record for entry:', /' ', A80,
     *         /' with length = ',I7, ', and Accession number = ', A6,
     1         /' The line that should have been the terminator record',
     2          ' is:', /' ', A80, // )
    5 FORMAT ( //' A DEFINITION line was not found for the entry:',
     *          /' ', A80, / )
C
C
      Do 100 i = 0, MxSymP, 1
         Raw(i) = 0
  100    continue
      NS = 0
      LsTot = 0
C
  200 READ( NFILE, 1, END = 500 )  LOCLN
C
      IF( LOCLN( 1:5 ) .EQ. 'LOCUS' )    THEN
         READ( LOCLN( 23:29 ), '(I7)' )  LENGTH    !   No. of bases in data
         Name( NS+1 ) = LOCLN( 13 : 22 )
         READ( NFILE, 1 )  DEFLN
         IF( DEFLN(1:10)   .EQ.   'DEFINITION' )    THEN
            SEQDEF = DEFLN( 12:80 )
            READ( NFILE, 1 )  ACC
            IF( ACC(1:3) .EQ. '   ' )    THEN   !  More Definition lines
               DO 230 NP = 70, 1, -1
  230             IF( SEQDEF(NP:NP) .NE. ' ' )    GO TO 240
  240          FIRST = NP + 2
               DO 250 NP = 94 - FIRST, 13, -1            !   94 = 13 + 80 + 1
  250             IF( ACC(NP:NP) .EQ. ' ' )    GO TO 260
  260          LC = NP - 1
               IF( LC .GT. 12 )    THEN
                  NDEND = FIRST + LC - 13
                  NA = 12
                  DO 270 ND = FIRST, NDEND, 1
                     NA = NA + 1
  270                SEQDEF( ND:ND ) = ACC( NA:NA )
               ELSE
               END IF
  280          READ( NFILE, 1 )  ACC
               IF( ACC(1:9) .NE. 'ACCESSION' )    GO TO 280
C
            ELSE IF( ACC(1:9) .EQ. 'ACCESSION' )    THEN
C
            ELSE                                          !   error
               WRITE( TERMNL, 3 )  LOCLN,  SEQDEF
               STOP ' GBacc'
            END IF
C
C ***   Process Accession Number(s)
C 
            ACNUM = ACC( 13:18 )
C
         ELSE
            WRITE( TERMNL, 5 )  LOCLN
            STOP ' GBdef'
         END IF
C
C ***       Search for ORIGIN line which marks the beginning of the
C ***       nucleotide sequence data
C
  350    READ( NFILE, 1 )  ORIGIN
         IF( ORIGIN( 1:6 ) .NE. 'ORIGIN' )    GO TO 350
C
C ***    Now read the nucleotide sequence and translate it into numeric form
C
         jblock = ( Length + 59 ) / 60
         jbgn = LsTot - 59
         Do 410 jb = 1, jblock, 1
            jbgn = jbgn + 60
            jend = jbgn + 59
            If( jend .gt. ( Length + LsTot ))   jend = Length + LsTot        
            READ( NFILE, 2 )  ( LSEQ( I ), I = 1, 60 )
            i = 0
            DO 400 j = jbgn, jend, 1
               i = i + 1
               SEQ( j ) = Lets( ICHAR( LSEQ( i ) ) )
               Raw( Seq(j) ) = Raw( Seq(j) ) + 1
  400       CONTINUE
  410    continue
      NS = NS + 1
      LSN( NS ) = Length
      LsTot = LsTot + Length
C
C ***   Find "//" line which marks the end of each data entry, then start
C ***   the next entry.
C
         READ( NFILE, 1 )   TERM
         IF( TERM(1:2)   .EQ.   '//' )    GoTo 200
C
C ***    ERROR in GenBank data file, correct and rerun program
C
         WRITE( TERMNL, 4 )  Name(NS), LENGTH, ACNUM, TERM
         CLOSE( UNIT = NFILE, STATUS = 'KEEP' )
         STOP ' GBtrm'
      ELSE
         GO TO 200            !   Skip header lines until first LOCUS line
      END IF
C
  500 RETURN
      END


 
C
C
      SUBROUTINE RdEMBL( TERMNL, NFILE, NSP, TotSqP, LsTot, NS, SEQ,
     *                   LSN, Name, MxSymP, Raw, Lets)
C
C ***   Subroutine RdEMBL reads an EMBL formatted ASCII sequence file and
C       returns the nucleic acid sequences coded in numeric form packed
C       contiguously in array SEQ.  The lengths of the sequences are
C       returned in array LSN.  The EMBL IDentification codes are returned
C       in the character array Name,
C
C       TERMNL is the FORTRAN unit number of the standard output file, where
C       error messages should be sent.  NFILE is the FORTRAN unit number for
C       the EMBL sequence data file.
C
C
      IMPLICIT NONE
C
C ***   Passed Variable declarations ( entry rdembl {read embl})
C
      INTEGER        TERMNL, NFILE, NSP, TotSqP, LsTot, NS, MxSymP
      INTEGER        SEQ( TotSqP ), Raw( 0:27 ), LSN( NSP ),
     *               Lets( 0:127 )
      CHARACTER*6    ACNUM
      CHARACTER*15   NAME( NSP )
      CHARACTER*80   SEQDEF
      Logical        Fertig
C
C ***   Local Variable declarations
C
      INTEGER   NAMEP, NPCHRS
C
      PARAMETER(  NAMEP = 10, NPCHRS = 20 )
C
      INTEGER        C1, C2, DELTA, EXTENT, I, L, LASTC, LIMIT, LS,
     *               START, FINISH, Length, jblock, jb, jbgn, jend, j
      CHARACTER*1    PUNCT( 0:NPCHRS ), LSeq( 60 )
      Character*15   SqName
      CHARACTER*31   WORD
      CHARACTER*80   IDLINE, ACLINE, DELINE, TERM
      LOGICAL        EMPTY
C
      INTEGER        LASTCH, LSEMBL                 !   function declaration
      EXTERNAL       GETWRD, LASTCH, LSEMBL
C
C
C ***   PUNCT contains the legal punctuation characters for VAX FORTRAN.
C ***   The first element, 0, is a sentinal needed by the binary search (POSN)
C ***   subroutine.  The punctuation characters are in ASCII collating
C ***   sequence order.
C
      DATA  PUNCT / ' ', ' ', ' ', '!', '"', '$', '%', '&', '''', '(',
     +              ')', '*', '+', ',', '-', '.', '/', ':', '<',  '=',
     +              '>'                                               /
C 
C
    1 FORMAT ( A80 )
    2 FORMAT ( 1 ( 4X, 6 ( 1X, 10(A1:) ) ) )
    3 FORMAT ( //' The expected ACcession line was not found for',
     *           ' entry:', /' ', A73 )
    4 FORMAT ( //' The expected DEfinition line was not found for',
     *           ' entry:', /' ', A73,
     1          /' with accession number:  ', A6 )
    5 FORMAT ( //' No TERMINATOR record "//" for entry: ', A10, I8, A6,
     *          /' ', A80, //' The line found instead of the expected',
     1           ' terminator was:', /' ', A72, // )
C
C
      PUNCT( 0 ) = CHAR( 0 )
      PUNCT( 1 ) = CHAR( 9 )
C
      Do 100 i = 0, MxSymP, 1
         Raw(i) = 0
  100    continue
      NS = 0
      LsTot = 0
C
  200 READ( NFILE, 1, END = 500 )  IDLINE
      FERTIG = .FALSE.
C
      IF( IDLINE( 1:2 )  .EQ.  'ID' )    THEN
         SqName(1:10) = IDLINE( 6:15 )
         SqName(11:15) = '     '
  220    READ( NFILE, 1 )  ACLINE
         IF( ACLINE( 1:2 )  .EQ.  'AC' )    THEN
         ELSE IF( ACLINE( 1:2 )  .EQ.  'XX'   .OR.
     *            ACLINE( 1:2 )  .EQ.  'CC' )    THEN
            GO TO 220
         ELSE 
            WRITE( TERMNL, 3 )  IDLINE( 1:73 )
            STOP 'IdLine'
         END IF
C         
         CALL GETWRD( 5, 73, ACLINE, START, FINISH, NPCHRS, EMPTY, WORD,
     *                PUNCT )
         ACNUM = WORD(1:6)
C
C ***   Skip the first line not of type "AC", it should be one of types
C ***   "DT", "XX", or "CC".  --  Read until line type "DEfinition" is found.
C
  290    READ( NFILE, 1 )  DELINE
         IF( DELINE( 1:2 )  .EQ.  'DE' )    THEN
            SEQDEF = DELINE( 6:72 )
         ELSE IF( DELINE( 1:2 )  .EQ.  'XX'   .OR.
     *            DELINE( 1:2 )  .EQ.  'CC'   .OR.
     1            DELINE( 1:2 )  .EQ.  'DT'   .OR.
     2            DELINE( 1:2 )  .EQ.  'AC' )    THEN
            GO TO 290
         ELSE
            WRITE( TERMNL, 4 )  IDLINE( 1:73 ), ACNUM
         END IF
C
         READ( NFILE, 1 )  DELINE
         IF( DELINE( 1:2 )  .EQ.  'DE' )    THEN
            LASTC = LASTCH( SEQDEF, 68 )
            DELTA = 79 - LASTC
            EXTENT = 6 + DELTA
            DO 310 I = EXTENT, 6, -1
               IF( DELINE(I:I)   .EQ.   ' '  )    THEN
                  LIMIT = I
                  GO TO 320
               ELSE
               END IF
  310          CONTINUE
            LIMIT = EXTENT - 1
  320       C1 = LASTC + 2
            C2 = C1 + LIMIT - 6
            SEQDEF( C1 : C2 ) = DELINE( 6 : LIMIT )
         ELSE
         END IF
C
C ***    Now find the SeQuence line and read the length of the sequence
C
  350    IF( DELINE(1:2) .NE. 'SQ' )    THEN
            READ( NFILE, 1 )  DELINE
            GO TO 350
         ELSE
         END IF
C
         LENGTH = LSEMBL( DELINE, 72 )
C
C ***   Now read the nucleotide sequence
C
         jblock = ( Length + 59 ) / 60
         jbgn = LsTot - 59
         Do 410 jb = 1, jblock, 1
            jbgn = jbgn + 60
            jend = jbgn + 59
            If( jend .gt. ( Length + LsTot ))    jend = Length + LsTot
            READ( NFILE, 2 )  ( LSEQ( I ), I = 1, 60 )
            i = 0
            DO 400 j = jbgn, jend, 1
               i = i + 1
               SEQ( j ) = Lets( ICHAR( LSEQ( i ) ) )
               Raw( Seq(j) ) = Raw( Seq(j) ) + 1
  400       CONTINUE
  410    continue
      NS = NS + 1
      LSN( NS ) = Length
      LsTot = LsTot + Length
      Name( NS ) = SqName
C
C ***   Now check to insure that the next line is a terminator line - then
C ***   go onto the next entry in the file.
C
         READ( NFILE, 1 )  TERM
         IF( TERM(1:2)   .EQ.   '//' )    GoTo 200
C
C ***   Error in EMBL file - correct and re-run the program
C
         WRITE( TERMNL, 5 )  SQNAME, LENGTH, ACNUM, SEQDEF, TERM(1:72)
         GO TO 500
C
      ELSE
         GO TO 200
      END IF
C
  500 FERTIG = .TRUE.
      RETURN
      END



C
C
      SUBROUTINE RdPIR( TERMNL, NFILE, NSP, TotSqP, LsTot, NS, SEQ,
     *                  LSN, Name, MxSymP, Raw, Lets)
C
C ***   Subroutine RdPIR reads an NBRF/PIR formatted ASCII sequence file and
C       returns the sequences, coded in numeric form, stored contiguously
C       in array SEQ.  The lengths of the sequences are returned in array
C       LSN.  The NBRF identification code is returned in the character
C       array Name.
C       NFILE is the FORTRAN unit number for the GenBank sequence data file.
C
      IMPLICIT NONE
C
C ***   Passed Variable declarations
C
      INTEGER        Termnl, NFILE, LsTot, TotSqP, NSP, NS, MxSymP
      INTEGER        SEQ( TotSqP ), Raw( 0:MxSymP ), LSN( NSP )
      INTEGER        Lets( 0:127 )
      CHARACTER*15   Name( NSP )
C
C ***   Local variable declarations
C
      INTEGER        BASE, L, LEN, LS, i
      CHARACTER*80   LocLn, DefLn
      CHARACTER*130  SqLine
C
      INTEGER        LASTCH                  !   function declaration
C
      EXTERNAL       LASTCH
C
C
    1 FORMAT ( A80 )
    2 FORMAT ( A130 )
C
C
      Do 80 i = 0, MxSymP, 1
         Raw(i) = 0
   80    continue
      NS = 0
      LsTot = 0
C
  100 READ( NFILE, 1, END = 500 )  LOCLN
C
      IF( LOCLN(1:1) .EQ. '>'   .AND.   LOCLN(4:4) .EQ. ';' )    THEN
         Name( NS+1 ) = LOCLN( 5:19 )
         READ( NFILE, 1 )  DEFLN
C
C ***    Now read the nucleotide sequence and translate it into numeric form
C
         LS = 0
  200    READ( NFILE, 2 )  SQLINE
         LEN = LASTCH( SQLINE, 130 )
         DO 300 L = 1, LEN
            BASE = Lets( ICHAR( SQLINE( L:L ) ) )
            IF((BASE .GE. 1 .AND. BASE .LE. 26) .or. BASE .eq. 28) THEN
               LS = LS + 1
               LsTot = LsTot + 1
               SEQ( LsTot ) = BASE
               Raw( BASE ) = Raw( BASE ) + 1
            ElseIf( Base .eq. 27 )    Then
               If( SqLine(L:L) .eq. '-' )    Then
                  LS = LS + 1
                  LsTot = LsTot + 1
                  SEQ( LsTot ) = BASE
                  Raw( BASE ) = Raw( BASE ) + 1
               ElseIf( SqLine(L:L) .eq. '.' )    Then
               EndIf
            ELSE IF( BASE .EQ. 100 )    THEN
               NS = NS + 1
               LSN( NS ) = LS
               GoTo 100
            ELSE
            END IF
  300       CONTINUE
         GO TO 200
C
      ELSE
         GO TO 100            !   Skip lines until first line of entry
      END IF
C
  500 RETURN
      END



C
C
      SUBROUTINE RFasta( TERMNL, SeqFl, NSP, TotSqP, LsTot, NS, SEQ,
     *                   LSN, Name, MxSymP, Raw, Lets)
C
C ***   Subroutine Fasta reads a file of sequences in Fasta format.
C ***   It places all of the sequences contiguously into a one
C ***   dimensional array after converting them to a numeric code.
C
      Implicit None
C
      Integer       Termnl, SeqFl, NSP, TotSqP, LsTot, NS, ls, MxSymP
      Integer       Seq( TotSqP ), LSN( NSP ), Raw( 0:MxSymP ),
     *              Lets( 0:127 )
      Character*15  Name( NSP )
C
      Integer        SeqLen, N, L, SeqChr
      Character*132  Line
C
      Integer        LastCh                         !   Function Declaration
C
C
    1 Format( A132 )
C
      NS = 0
      lstot = 0
      Do 50 N = 0, MxSymP, 1
         Raw(n) = 0
   50    continue
C
C ***   Read the sequence file and count all of the alphabetic characters
C ***   on all lines that are not a "beginning of sequence" record.
C
  100 Read( SeqFl, 1, End = 500 )  Line
      If( Line(1:1) .eq. '>' )    Then
         If( NS .gt. 0 )    LSN( NS ) = ls
         NS = NS + 1
         Name(NS) = line(2:16)
         ls = 0
      Else
         SeqLen = Lastch( Line, 132 )
         If( SeqLen .eq. 0 )    GoTo 100
         Do 150 L = 1, SeqLen
            SeqChr = Lets( Ichar( Line(l:l) ) )
            If( SeqChr .ge. 1   .and.   SeqChr .le. MxSymP )    Then
               lstot = lstot + 1
               Seq( lstot ) = SeqChr
               Raw( SeqChr ) = Raw( SeqChr ) + 1
               ls = ls + 1
            Else
            EndIf
  150       continue
      EndIf
      GoTo 100
C
  500 LSN( NS ) = ls
C
      Return
      End



C
C
      SUBROUTINE MolGn( TERMNL, NFILE, NSP, TotSqP, LsTot, NS, SEQ,
     *                  LSN, Name, MxSymP, Raw, Lets)
C
C ***   Subroutine MOLGN reads a MOLGEN/STANFORD formatted ASCII sequence
C       file and returns the nucleic acid sequences, coded in numeric form,
C       stored contiguously in array SEQ.  The lengths of the sequences
C       are returned in Array LSN.  The MOLGEN identification codes are
C       returned in the character array Names.
C       NFILE is the FORTRAN unit number for the GenBank sequence data file.
C
      IMPLICIT NONE
C
C ***   Passed Variable declarations
C
      INTEGER        Termnl, NFILE, NSP, TotSqP, LsTot, NS, MxSymP
      INTEGER        SEQ( TotSqP ), Raw( 0:MxSymP ), Lets( 0:127 ),
     *               LSN( NSP )
      CHARACTER*15   Name( NSP )
      CHARACTER*80   DEFLN
C
C ***   Local variable declarations
C
      INTEGER        BASE, L, LEN, LS, i
      INTEGER        NUC( 0:127 )
      CHARACTER*80   LOCLN
      CHARACTER*130  SQLINE
C
      INTEGER        LASTCH                  !   function declaration
C
      EXTERNAL       LASTCH
C
C
    1 FORMAT ( A80 )
    2 FORMAT ( A130 )
C 
C
      Do 80 i = 0, MxSymP, 1
         Raw(i) = 0
   80    continue
      NS = 0
      LsTot = 0
C
  100 READ( NFILE, 1, END = 500 )  DEFLN
C
      IF( DEFLN(1:1)   .EQ.   ';' )    THEN
         DEFLN(1:1) = ' '
  110    READ( NFILE, 1, END = 500 )  LOCLN
         IF( LOCLN(1:1) .EQ. ';' )    GO TO 110            !   skip comments
         Name( NS+1 ) = LOCLN( 1:15 )
C
C ***    Now read the nucleotide sequence and translate it into numeric form
C
         LS = 0
  200    READ( NFILE, 2 )  SQLINE
         LEN = LASTCH( SQLINE, 130 )
         DO 300 L = 1, LEN
            BASE = Lets( ICHAR( SQLINE( L:L ) ) )
            IF( BASE .GE. 1   .AND.   BASE .LE. MxSymP )   THEN
               LS = LS + 1
               LsTot = LsTot + 1
               SEQ( LsTot ) = BASE
               Raw( Seq(LsTot) ) = Raw( Seq(LsTot) ) + 1
            ELSE IF( BASE .EQ. 201   .or.   Base .eq. 202 )    THEN
               NS = NS + 1
               LSN( NS ) = LS
               GoTo 100
            ELSE
            END IF
  300       CONTINUE
         GO TO 200
C
      ELSE
         GO TO 100            !   Skip lines until first line of entry
      END IF
C
  500 RETURN
      END



C
C
      SUBROUTINE Simpl( TERMNL, NFILE, NSP, TotSqP, LsTot, NS, SEQ,
     *                  LSN, Name, MxSymP, Raw, Lets)
C
C ***   Subroutine SIMPL reads a simple ASCII sequence file and returns
C       the sequences, coded in numeric form, stored contiguously in
C       array SEQ.  The length of the sequences are returned in array
C       LSN.  The first 10 characters of the identification lines are
C       returned in the character array Name.
C       NFILE is the FORTRAN unit number for the GenBank sequence data file.
C
      IMPLICIT NONE
C
C ***   Passed Variable declarations
C
      INTEGER        Termnl, NFILE, NSP, TotSqP, LsTot, MxSymP, NS
      INTEGER        SEQ( TotSqP ), LSN( NSP ), Raw( 0:MxSymP ),
     *               Lets( 0:127 )
      CHARACTER*15   Name( NSP )
C
C ***   Local variable declarations
C
      INTEGER        BASE, L, LEN, LS, i
      INTEGER        NUC( 0:127 )
      CHARACTER*80   DEFLN
      CHARACTER*130  SQLINE
C
      INTEGER        LASTCH                  !   function declaration
C
      EXTERNAL       LASTCH
C
C
    1 FORMAT ( A80 )
    2 FORMAT ( A130 )
C 
C
      Do 80 i = 0, MxSymP, 1
         Raw(i) = 0
   80    continue
      NS = 0
      LsTot = 0
C
  100 READ( NFILE, 1, END = 500 )  DEFLN
      LEN = LASTCH( DEFLN, 80 )
      IF( LEN .EQ. 0 )    GO TO 100
      Name( NS+1 ) = DEFLN( 1:15 )
C
C ***    Now read the nucleotide sequence and translate it into numeric form
C
      LS = 0
  200 READ( NFILE, 2, END = 400 )  SQLINE
      LEN = LASTCH( SQLINE, 130 )
      IF( LEN .EQ. 0 )    THEN
         NS = NS + 1
         LSN( NS ) = LS
         GoTo 100
      ELSE
         DO 300 L = 1, LEN
            BASE = Lets( ICHAR( SQLINE( L:L ) ) )
            IF( BASE .GE. 1   .AND.   BASE .LE. MxSymP )   THEN
               LS = LS + 1
               LsTot = LsTot + 1
               SEQ( LsTot ) = BASE
               Raw( Base ) = Raw( Base ) + 1
            ELSE
            END IF
  300       CONTINUE
         GO TO 200
      END IF
C
  400 NS = NS + 1
      LSN( NS ) = LS
  500 RETURN
      END



C
C
      Subroutine FixSeq( LsTot, TotSqP, Seq, MxSymP, Raw, Comp, Gap,
     *                   MaxLet, ASize, SqPrnt, n25, ModAA, Protin,
     1                   DNA, BackFr, UnKnown )
C
C ***   Subroutine FixSeq determines whether the sequences are protein
C ***   or nucleic acid and translates the original numeric code based
C ***   on the complete alphabet based on either the protein or nucleic
C ***   acid alphabet.
C
      Implicit None
C
      Integer       LsTot, TotSqP, MaxLet, Gap, MxSymP, n25, ModAA,
     *              ASize, UnKnown
      Integer       Seq( TotSqP ), Raw( 0:MxSymP ), Comp( 0:MxSymP )
      Real*4        BackFr( MxSymP )
      Character*1   SqPrnt(n25)
      Logical       Protin, DNA
C
      Integer        N, PrOnly, ACGT, Total
      Integer        NucMap(27), AAMap(27)
      Real*4         FTotal
      Character*1    AA1(24), Base(6)
C
C
      Data  Base   / 'a', 'c', 'g', 't', 'n', '-' /
C
      Data  NucMap /  1,  0,  2,  0,  0,  0,  3,  0,  0,  0,  0,  0,
     *                0,  5,  0,  0,  0,  5,  0,  4,  4,  0,  0,  5,
     1                5,  0,  6  /
C
C      Data  AA3  / 'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu',
C     *             'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe',
C     1             'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 'Asx',
C     2             'Glx', 'Unk'  '---'  /
C
C
      Data  AA1   / 'a', 'r', 'n', 'd', 'c', 'q', 'e', 'g', 'h', 'i',
     *              'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v',
     1              'b', 'z', 'x', '-'   /
C
      Data  AAMap  /  1, 21,  5,  4,  7, 14,  8,  9, 10,  0, 12, 11,
     *               13,  3,  0, 15,  6,  2, 16, 17,  0, 20, 18, 23,
     1               19, 22, 24  /
C
C
C ***   Determine whether the sequences are Nucleic Acid or Protein by
C ***   counting those letters that code for amino acids but not any of
C ***   the IUB nucleotides ( E, F, I, L, P, Q, and Z )
C
  500 PrOnly = Raw(5) + Raw(6) + Raw(9) + Raw(12) + Raw(16) + Raw(17)
     *       + Raw(26)
      ACGT = Raw(1) + Raw(3) + Raw(7) + Raw(20) + Raw(21)
C
      If( PrOnly .gt. 0 )    Then
         Protin = .True.
         DNA = .False.
         MaxLet = 23
         ASize = 20
         Gap = 24
         UnKnown = 23
         Comp(0) = Raw(23) + Raw(22) + Raw(2)
         Comp(1) = Raw(1)
         Comp(2) = Raw(18)
         Comp(3) = Raw(14)
         Comp(4) = Raw(4)
         Comp(5) = Raw(3)
         Comp(6) = Raw(17)
         Comp(7) = Raw(5)
         Comp(8) = Raw(7)
         Comp(9) = Raw(8)
         Comp(10) = Raw(9)
         Comp(11) = Raw(12)
         Comp(12) = Raw(11)
         Comp(13) = Raw(13)
         Comp(14) = Raw(6)
         Comp(15) = Raw(16)
         Comp(16) = Raw(19)
         Comp(17) = Raw(20)
         Comp(18) = Raw(23)
         Comp(19) = Raw(25)
         Comp(20) = Raw(22)
         Comp(28) = Raw(28)
         Comp(ModAA) = Comp( ModAA ) + Raw( 28 )
         Do 540 n = 1, LsTot
            Seq(n) = AAMap( Seq(n) )
  540       continue
         Do 550 n = 1, Gap, 1
            SqPrnt(n) = AA1(n)
  550       continue
C
      ElseIf( ACGT .gt. 0 )    Then
         Protin = .False.
         MaxLet = 5
         ASize = 4
         Gap = 6
         UnKnown = 5
         Comp(0) = Raw(14) + Raw(24)
         Comp(1) = Raw(1)
         Comp(2) = Raw(3)
         Comp(3) = Raw(7)
         Comp(4) = Raw(20) + Raw(21)
         Comp(28) = Raw(28)
         Comp(ModAA) = Comp( ModAA ) + Raw( 28 )
         If( Raw(20) .ge. Raw(21) )    Then
            DNA = .True.
         Else
            DNA = .False.
            Base(4) = 'u'
         EndIf
         Do 560 n = 1, LsTot
            Seq(n) = NucMap( Seq(n) )
  560       continue
         Do 570 n = 1, Gap, 1
            SqPrnt(n) = Base(n)
  570       continue
      EndIf
C
      Total = 0
      Do 700 n = 1, ASize, 1
         Total = Total + Comp(n)
         BackFr(n) = Float( Comp(n) )
  700    continue
      FTotal = Float( Total )
      Do 750 n = 1, ASize
         BackFr(n) = BackFr(n) / FTotal
  750    continue
C
      Return
      End



C
      Subroutine GetLst( Termnl, KeyBrd, Count, NameP, NPChrs, MaxLin,
     *                   ILine, Names, Punct )
C
C ***   Subroutine GetLst gets a list of words from a line of input, ILINE,
C ***   and stores the individual words in an array NAMES.  If the last
C ***   character in the line is < another line of is solicited to include
C ***   in the list (NAMES).  A word is defined as a string of characters
C ***   delimited by the characters in the punctuation list (PUNCT).  The
C ***   maximum number of characters allowed on a line is set by MAXLIN.
C ***   COUNT returns the number of words found.
C
      IMPLICIT NONE
C
      INTEGER        TERMNL, KEYBRD, COUNT, MAXLIN, NAMEP, NPCHRS
      CHARACTER*1    PUNCT( 0 : NPCHRS )
      CHARACTER*(*)  NAMES( NAMEP )
      CHARACTER*(*)  ILINE
C
      INTEGER       I, K, IBGN, IEND, NEXTC, LASTC
C      CHARACTER*10  FMT2
      CHARACTER*31  WORD
      LOGICAL       EMPTY
C
      INTEGER   LASTCH                                    !   function name
      EXTERNAL  GETWRD
C      EXTERNAL  AFIELD
C
C
C      DATA  FMT2  /  '(A       )'  /            !   3 : 9 out of 10
C
C
    1 FORMAT ( //' Enter an additional line to be processed.  Enter',
     *         ' an < at the end of',
     1        /' the line if you need to continue the input onto',
     2         ' another line.',//)
C
C    3 FORMAT ( //'     **********     WARNING     **********',
C     *          /' A completely blank line was entered - please',
C     1           ' respond with valid input.'
C     2          /' A second blank line will cause the program to',
C     3           ' ABORT.'//)
C
C
C      CALL AFIELD( TERMNL, MAXLIN, 7, FMT2(3:9) )
      COUNT = 0
      NEXTC = 1
      LASTC = LASTCH( ILINE, MAXLIN )
C
  100 CALL GETWRD( NEXTC, LASTC, ILINE, IBGN, IEND, NPCHRS, EMPTY,
     *             WORD, PUNCT )
C
      IF ( .NOT. EMPTY )    THEN
C
         COUNT = COUNT + 1
         NAMES( COUNT ) = WORD
         NEXTC = IEND + 1
         IF( IEND .LT. LASTC )    GO TO 100
C
C      ELSE IF( EMPTY   .AND.   COUNT .GT. 0 )    THEN
C
C         IF( ILINE(LASTC:LASTC) .EQ. '<' )    THEN
C            WRITE( TERMNL, 1 )
C            READ( KEYBRD, FMT2 )  ILINE
C            NEXTC = 1
C            LASTC = LASTCH( ILINE, MAXLIN )
C            GO TO 100
C         ELSE
C         END IF
C
      ELSE IF( EMPTY   .AND.   COUNT .LE. 0 )    THEN
C
C ***   The following code which has been inactivated by making it comment
C ***   statements is suitable for interactive use.  It traps and resolicits
C ***   blank input lines.  If you decide to use this code FORMAT 3 must
C ***   also be activated by removing the C in position 1 of each line.
C
C         WRITE( KEYBRD, 3 )
C         READ( KEYBRD, FMT2 )  ILINE
C         LASTC = LASTCH( ILINE, MAXLIN )
C         IF( LASTC .EQ. 0 )   STOP ' invalid input.'
C         NEXTC = 1
C         GO TO 100
C
         RETURN
C
      END IF
C
      RETURN
      END
C
C
C
      Subroutine GetWrd( BGN, LAST, LINE, START, FINISH, NPCHRS, EMPTY,
     *                   WORD, PUNCT )
C
C ***   Subroutine GetWrd gets the next FORTRAN symbolic name from a line
C ***   of text.  The subroutine begins searching at position BGN of the
C ***   line of text and stops at position LAST.  If no symbolic name is
C ***   found in this region of the line the logical variable EMPTY is set
C ***   .TRUE..  If a variable is found EMPTY is set .FALSE. and the symbolic
C ***   name is transfered to the variable WORD for return to the calling
C ***   program.  If a symbolic name is found, the variables START and FINISH
C ***   mark positions of the first and last characters of the name on the
C ***   line of text (LINE).  Words longer than 31 characters are truncated
C ***   to 31 characters - but this limit can be changed by redeclaring
C ***   WORD and resetting the LIMIT parameter.
C
      IMPLICIT NONE
C
      INTEGER        LIMIT, LIMITM
C
      PARAMETER  (   LIMIT = 31,   LIMITM = LIMIT - 1 )
C
      INTEGER        BGN, LAST, START, FINISH, NPCHRS
      CHARACTER*1    PUNCT( 0:NPCHRS )
      CHARACTER*31   WORD
      CHARACTER*(*)  LINE
      LOGICAL        EMPTY
C
      INTEGER   I, LENGTH
      INTEGER   POSN                                       !   function name
C
      EXTERNAL  POSN
C
C
      START = 0
      FINISH = 0
      WORD = '                               '
C
      IF( BGN .GT. LAST )    THEN
         EMPTY = .TRUE.
         RETURN
      ELSE
      END IF
C
C ***   Find the first non-punctuation character on the line
C
      DO 120 I = BGN, LAST
         IF( POSN( NPCHRS, LINE(I:I), PUNCT ) .EQ. 0 )    THEN
            START = I
            EMPTY = .FALSE.
            GO TO 150
         ELSE
         END IF
  120    CONTINUE
      EMPTY = .TRUE.
      RETURN
C
C ***   Find the first punctuation character after the START location.
C
  150 DO 200 I = START + 1, LAST
         IF( POSN( NPCHRS, LINE(I:I), PUNCT ) .NE. 0 )    THEN
            FINISH = I - 1
            LENGTH = FINISH - START + 1
            WORD(1:LENGTH) = LINE(START:FINISH)
            RETURN
         ELSE
         END IF
  200    CONTINUE
C
      IF( LAST - START  .LE. LIMITM )    THEN
         FINISH = LAST
      ELSE
         FINISH = START + LIMITM
      END IF
      LENGTH = FINISH - START + 1
      WORD(1:LENGTH) = LINE(START:FINISH)
C
      RETURN
      END
C
C
      INTEGER FUNCTION LSEMBL( TEXT, LENGTH )
C
C  ****************************************************************************
C  **      Copyright Pittsburgh Supercomputing Center  --  June 1988         **
C  ***************************************************************************
C
C ***   LSEMBL returns the numeric value of the first string of contiguous
C ***   digits beginning after position 12 in the character variable TEXT.
C ***   The function returns only positive values since negative values are
C ***   not expected (the number being translated is the count of bases in
C ***   the sequence of an entry in the EMBL data file format).  Zero is
C ***   returned if no digits are encountered on the line before either
C ***   the end of the line or the string BP is encountered.
C
      IMPLICIT  NONE
C
      INTEGER        ZEROCH
C
      PARAMETER  (  ZEROCH = 48  )      !   location of zero "0" in ASCII
C                                           collating sequence
C
      INTEGER        LENGTH
      CHARACTER*(*)  TEXT
C
      INTEGER        L, FIRST, NUMBER
C
      DO 100 L = 13, LENGTH
         IF( TEXT( L:L ) .EQ. ' ' )    THEN
         ELSE IF( LGE( TEXT(L:L), '0') .AND. LLE( TEXT(L:L), '9') ) THEN
            FIRST = L
            NUMBER = ICHAR( TEXT( L:L ) ) - ZEROCH
            GO TO 150
         ELSE IF(TEXT(L:L+1) .EQ. 'BP' .OR. TEXT(L:L+1) .EQ. 'bp')  THEN
            LSEMBL = 0
            RETURN
         ELSE
         END IF
  100    CONTINUE
C
      LSEMBL = 0
      RETURN
C
  150 DO 200 L = FIRST + 1, LENGTH
         IF( LGE( TEXT(L:L), '0' )  .AND.  LLE( TEXT(L:L), '9' ) )  THEN
            NUMBER = ( 10 * NUMBER ) + ( ICHAR( TEXT(L:L) ) - ZEROCH )
         ELSE
            LSEMBL = NUMBER
            RETURN
         END IF
  200    CONTINUE
C
      LSEMBL = NUMBER
      RETURN
      END
C
C
      SUBROUTINE WHATFM( SEQUNT, TERMNL, FLFMT )
C
C ***   Subroutine WHATFM determines whether a sequence file is EMBL format
C ***   GenBank format, Genbank format without a header, or an undetermined
C ***   format.  The unit number of the file to be checked is passed to the
C ***   subroutine in SEQUNT.  TERMNL is the unit number to which normal
C ***   results (and error messages) should be sent.  The file format is
C ***   returned as an integer value in the variable FLFMT.
C
C         FLFMT = 0   Undetermined file format, not GenBank or EMBL
C         FLFMT = 1   GenBank file format
C         FLFMT = 2   EMBL file format
C         FLFMT = 3   NBRF file format
C         FLFMT = 4   MolGen/Stanford file format
C         FLFMT = 5   Fasta (Pearson) format
C
      IMPLICIT NONE
C
C ***   Passed variables
C
      INTEGER      SEQUNT, TERMNL, FLFMT
C
C ***   Local Variables and constants
C
      INTEGER      GENBNK, EMBL, NBRF, FASTA, MOLGEN, Simple, NONE,
     *             GCGMSF
C
      PARAMETER ( NONE = 0, GENBNK = 1, EMBL = 2, NBRF = 3, MOLGEN = 4,
     *            FASTA = 5, Simple = 6, GCGMSF = 7 )
C
C
      INTEGER       Len, I, HEADER, NBlank
      CHARACTER*5   LOCUS
      CHARACTER*12  DBANK( 0:7 )
      CHARACTER*70  GBHEAD( 9 )
      CHARACTER*80  LINE, TEST
C
      Integer       LastCh
C
      DATA  DBANK  /  'Undetermined', 'GenBank     ', 'EMBL        ',
     *                'NBRF        ', 'MolGen      ', 'FASTA       ',
     1                'Simple      ', 'GCG MSF     ' /
C
      EXTERNAL   ALLCAP
C
    1 FORMAT ( A80 )
    2 FORMAT ( //' The format of sequence data file is:  ', A12 )
    3 FORMAT ( /'                   GenBank file header information',/)
    4 FORMAT (  ' ', A70 )
    5 FORMAT ( /' The header information was not present in the',
     *          ' GenBank sequence file.' )
C
C
      READ( SEQUNT, 1 )  LINE
      TEST = LINE
      CALL ALLCAP( TEST, 80 )
      HEADER = 0
      NBlank = 0
C 
C
      IF( TEST(1:2)   .EQ.   'ID' )    THEN
         FLFMT = EMBL
      ELSE IF( TEST(1:5)   .EQ.   'LOCUS' )    THEN
         FLFMT = GENBNK
      ELSE IF( TEST(20:45)  .EQ.  'GENETIC SEQUENCE DATA BANK' )   THEN
         FLFMT = GENBNK
  200    HEADER = HEADER + 1
         GBHEAD( HEADER ) = LINE( 1:70 )
         READ( SEQUNT, 1 )  LINE
         CALL ALLCAP( LINE, 5 )
         IF( LINE(1:5) .NE. 'LOCUS'  .AND.  HEADER .LT. 9 )   GO TO 200
      ELSE IF( TEST(1:1) .EQ. '>' )    THEN
         IF( TEST(4:4) .EQ. ';' )    THEN
            FLFMT = NBRF
         ELSE
            FLFMT = FASTA
         END IF
      ELSE IF( TEST(1:1) .EQ. ';' )    THEN
         FLFMT = MOLGEN
      ELSE
         Do 300 i = 1, 2500
         If( Index( Line, 'MSF: ' )   .gt. 0    .and.
     *       Index( Line, 'Check: ' ) .gt. 0    .and.
     *       Index( Line, '..' ) .gt.  0  )         Then
            FlFmt = GCGMSF
            GoTo 500
         ElseIf( LastCh( Line, 80 ) .eq. 0 )    Then
            NBlank = NBlank + 1
         Else
         EndIf
         READ( SeqUnt, 1, End = 310 )  Line
  300       continue
  310    If( NBlank .gt. 0 )    Then
            FlFmt = Simple
         Else
            FLFMT = None
         EndIf
      END IF
C
  500 REWIND( SEQUNT )
      WRITE( TERMNL, 2 )  DBANK( FLFMT )
      IF( FLFMT .EQ. GENBNK   .AND.   HEADER .GT. 0 )    THEN
         WRITE( TERMNL, 3 )
         WRITE( TERMNL, 4 )  ( GBHEAD( I ), I = 1, HEADER )
      ELSE IF( FLFMT .EQ. GENBNK   .AND.   HEADER .EQ. 0 )    THEN
         WRITE( TERMNL, 5 )
      ELSE
      END IF
C
      RETURN
      END
C
C
      Subroutine GetQij( Termnl, KeyBrd, OutFl, Bits, Qij, QSum, n20,
     *                   n24, MatAve, MaxMat, MinMat, MatSiz, MatDef,
     1                   SqPrnt )
C
      Implicit None
C
      Integer       KeyBrd, Termnl, OutFl, i, j, n20, n24, PamID,
     *              MatSiz
      Real*4        Qij(n20,n20), QSum(n20), MatAve, MaxMat, MinMat
      Character*1   SqPrnt(n24)
      Character*4   Bits
      Character*60  MatDef
C
      Character  AAOrd(20)
c      Integer  AAOrd(20)
      Data  AAOrd / 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
     *              'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' /
C
C
    1 Format(//' Please select, by number, one of the amino acid',
     *         ' transition frequency',
     1        /' matrices "q(i,j)" from the menu below:',
     2  //'   Amino Acid transition, Blosum percent clustering.',
     3   /'       1.  Blosum_30      2.  Blosum_40      3.  Blosum_50',
     4   /'       4.  Blosum_60      5.  Blosum_70      6.  Blosum_80',
     5   /'       7.  Blosum_90      8.  Blosum_100',
     6  //'   Amino Acid transition, PAM - accepted point mutations.',
     7   /'       9.  PAM_80        10.  PAM_120       11.  PAM_160',
     8   /'      12.  PAM_200       13.  PAM_250       14.  PAM_300',
     9  //'   DNA similarity, equal transitions and transversions.',
     A   /'      15.  PAM_20        16.  PAM_50        17.  PAM_80',
     B   /'      18.  PAM_110',
     C  //'   DNA similarity, 3 to 1 transitions to transversions.',
     D   /'      19.  PAM_20        20.  PAM_50        21.  PAM_80',
     E   /'      22.  PAM_110',
     F  //'   23.  Read a custom, transition frequency matrix from a',
     G    ' file.', // )
    2 Format(/' Qij matrix -  Maximum =', F7.4, ',   Minimum =', F7.4,
     *        ',   Average=', F7.4 )
    3 Format(' QSum( ',A1,' ) = ', F8.4 )
    6 Format( ' ', 20F4.3: )
    7 Format( ' ', 20( '  ',A1,' ') )
    8 Format(/' Bits:(',A4,')',/' MatDef:', A60 )
C
C
  100 Write( Termnl, 1 )
      Read( KeyBrd, * )   PamID
C
      If( PamId .eq. 1 )    Then
         Call Bl30( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 2 )    Then
         Call BL40( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 3 )    Then
         Call BL50( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 4 )    Then
         Call BL60( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 5 )    Then
         Call BL70( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 6 )    Then
         Call BL80( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 7 )    Then
         Call BL90( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 8 )    Then
         Call BL100( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 9 )    Then
         Call Pam80( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 10 )    Then
         Call Pam120( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 11 )    Then
         Call Pam160( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 12 )    Then
         Call Pam200( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 13 )    Then
         Call Pam250( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 14 )    Then
         Call Pam300( Bits, Qij, n20, MatDef )
         MatSiz = 20
      ElseIf( PamID .eq. 15 )    Then
         Call dp20( Bits, Qij, n20, MatDef )
         MatSiz = 4
      ElseIf( PamID .eq. 16 )    Then
         Call dp50( Bits, Qij, n20, MatDef )
         MatSiz = 4
      ElseIf( PamID .eq. 17 )    Then
         Call dp80( Bits, Qij, n20, MatDef )
         MatSiz = 4
      ElseIf( PamID .eq. 18 )    Then
         Call dp110( Bits, Qij, n20, MatDef )
         MatSiz = 4
      ElseIf( PamID .eq. 19 )    Then
         Call dp20b( Bits, Qij, n20, MatDef )
         MatSiz = 4
      ElseIf( PamID .eq. 20 )    Then
         Call dp50b( Bits, Qij, n20, MatDef )
         MatSiz = 4
      ElseIf( PamID .eq. 21 )    Then
         Call dp80b( Bits, Qij, n20, MatDef )
         MatSiz = 4
      ElseIf( PamID .eq. 22 )    Then
         Call dp110b( Bits, Qij, n20, MatDef )
         MatSiz = 4
      ElseIf( PamID .eq. 23 )    Then
         Call RdQij( Termnl, KeyBrd, Bits, Qij, n20, MatDef )
         MatSiz = 4
      Else
         GoTo 100
      EndIF
C
      MatAve = 0.0
      MinMat = 10000.0
      MaxMat = -10000.0
      Do 500 i = 1, MatSiz, 1
         QSum(i) = 0.0
         Do 490 j = 1, MatSiz, 1
            QSum(i) = QSum(i) + Qij(i,j)
            MatAve = MatAve + Qij(j,i)
            If( Qij(i,j) .gt. MaxMat )    Then
               MaxMat = Qij(i,j)
            ElseIf( Qij(i,j) .lt. MinMat )    Then
               MinMat = Qij(i,j)
            Else
            EndIf
  490       continue
  500    continue
      MatAve = MatAve / ( MatSiz * MatSiz )
      Write( Termnl, 2 )  MaxMat, MinMat, MatAve
      Write( OutFl, 2 )   MaxMat, MinMat, MatAve
      Write( OutFl, 8 )  Bits, MatDef
      Write( OutFl, 7 )  ( SqPrnt(i), i = 1, MatSiz )
      Do 600 i = 1, MatSiz
         Write( OutFl, 6 )  ( Qij(j,i), j = 1, MatSiz )
  600    continue
      Write( OutFl, 7 )  ( SqPrnt(i), i = 1, MatSiz )
      Do 650 i = 1, MatSiz, 1
         Write( OutFl, 3 )  SqPrnt(i), QSum(i)
  650    continue
C
      Return
C      Stop ' Transition Frequency Matrix defined.'
      End
C
C
      Subroutine RdQij( Termnl, KeyBrd, Bits, Qij, n20, MatDef )
C
C ***   Subroutine RdQij reads a file containing a transition frequency matrix
C ***   and reorders it to match the order of amino acids or nucleotides
C ***   used in the program.  The initial order is specified by the order
C ***   the characters in array InOrd (Input Order) which is the third
C ***   of the input file.  Either AAOrd or NucOrd depending on whether
C ***   the matrix is for Amino Acid or Nucleic Acid sequences specifies
C ***   the final order of the similarity matrix as it will be stored as
C ***   a symmetric matrix in array Qij.
C ***   The input file must have the following information and format.
C
C       *****     Input File Specifications     *****
C
C Line 1.  A title for the matrix, up to 60 characters of identifying
C          information which will be used to label the results.
C
C Line 2.  Up to four characters or digits that describe the scale or
C          entropy content of the frequencies.  This line can be blank
C          but it must be present.  As an example the Blosum matrices
C          for amino acid similarities are frequently on a scale of
C          1/3 bit of information per unit of score, thus 1/3 is an
C          appropriate notation in this case.
C
C Line 3.  The order of the transition frequency matrix that follows
C          in the file.  The order is specified by listing the single
C          letter codes of either the amino acids or nucleotides in
C          the order that they would be used to label the rows and
C          columns of the matrix.  For example, the amino acid order
C          used internally by the GroupEntropy program is based on the
C          alphabetical order of the three letter codes for amino acids.
C          Thus the line would be:
C
C          A R N D C Q E G H I L K M F P S T W Y V
C
C          The letters may be either upper or lower case and need not
C          be separated by spaces or other punctuation.
C
C Line 4+  The fourth and each subsequent line contains a complete row
C          of the matrix.  Thus each line will contains four nucleotide
C          or 20 amino acid transition frequencies.  The numbers must be
C          separated by spaces.  Each column of the matrix must contain
C          the data for the transitions from a single kind sequence
C          residue (i.e., a single amino acid such as alanine) to the
C          other possible kinds of sequence residue (i.e., the other
C          19 amino acids) and frequency with which it is not observed
C          to change.
C          The scores must be integers that are 10,000 times the actual
C          fractional value of the transition frequencies.
C
C
      Integer       Termnl, KeyBrd, i, j, MatFl, MSize, n20
      Integer       InMat(20,20), Map(20)
      Real*4        Qij(n20,n20)
      Character*1   Next, AAOrd(20), InOrd(20), NucOrd(4)
      Character*4   Bits
      Character*60  MatDef
      Character*80  MtFile, Line
      Logical       empty
C
      Data  AAOrd / 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
     *              'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' /
C
      Data  NucOrd  /  'A', 'C', 'G', 'T'  /
C
    1 Format(//' Please enter the name of the file holding the',
     *         ' complete',
     1        /' (square) transition frequency matrix.'/)
    2 Format( A80 )
    3 Format(//' Your transition frequency matrix is the wrong size.',
     *        /' Matrices for proteins must have 20 rows and columns',
     1        /' while nucleic acid matrices must have 4.  Yours',
     2        /' has ', I4, ' rows and columns.'//)
C    4 Format(/' The original order was:', / ' ', 20( A1,', ' ), / )
C    5 Format(/' The reordered matrix:', //' ', 20( '  ',A1,' '),/)
C
C
      MatFl = 2
      Write( Termnl, 1 )
      Read( KeyBrd, 2 )  MtFile
      Open( Unit = MatFl, file = MtFile, status = 'old' )
      Read( MatFl, 2 )   Line
      Call PakNam( Line, 80, empty )
      MatDef = Line(1:60)
      Read( MatFl, 2 )   Line
      Call PakNam( Line, 80, empty )
      Bits = Line(1:4)
      Read( MatFl, 2 )   Line
      MSize = 0
      Do 150 i = 1, 80
         If( Line(i:i) .ge. 'A' .and. Line(i:i) .le. 'Z' )    Then
            MSize = MSize + 1
            InOrd(MSize) = Line(i:i)
         ElseIf( Line(i:i) .ge. 'a' .and. Line(i:i) .le. 'z' )   Then
            MSize = MSize + 1
            InOrd(MSize) = Char( IChar( Line(i:i) ) - 32 )
         Else
         EndIf
  150    continue
      Do 200 i = 1, MSize, 1
         Read( MatFl, * )  ( InMat(j,i), j = 1, 20, 1 )
  200    continue
      Close( Unit = MatFl, status = 'keep' )
C
C ***   The matrix has been read and returned to a square shape,
C ***   now create a map of the needed reordering.
C
      If( MSize .eq. 20 )    Then
         Do 250 i = 1, 20
            Next = AAOrd(i)
            Do 240 j = 1, 20
               If( InOrd(j) .eq. Next )    Then
                  Map(i) = j
                  GoTo 250
               Else
               EndIf
  240          continue
  250       continue
      ElseIf( MSize .eq. 4 )    Then
         Do 280 i = 1, 4, 1
            If( InOrd(i) .eq. 'U' )    Then
               InOrd(i) = 'T'
            Else
            EndIf
  280       continue
         Do 300 i = 1, 4, 1
            Next = NucOrd(i)
            Do 290 j = 1, 4, 1
               If( InOrd(j) .eq. Next )    Then
                  Map(i) = j
                  GoTo 300
               Else
               EndIf
  290          continue
  300       continue
      Else
         Write( Termnl, 3 )  MSize
         Stop ' Illegal Matrix'
      EndIf
C
C ***   The map has been created, now reorder the matrix
C
      Do 500 i = 1, MSize
         Do 490 j = 1, MSize
            Qij(j,i) = InMat( Map(j), Map(i) )
  490       continue
  500    continue
C
C      Open( Unit = 3, file = 'Matrix.new', status = 'new',
C     *      access = 'Sequential', Form = 'formatted',
C     *      CarriageControl = 'list' )
C      Write( 3, 4 )  ( InOrd(i), i = 1, 20 )
C      Write( 3, 5 )  ( AAOrd(i), i = 1, 20 )
C      Do 600 i = 1, 23
C         Write( 3, 6 )  ( Qij(j,i), j = 1, 20 )
C  600    continue
C      Write( 3, 7 )  ( AAOrd(i), i = 1, 20 )
C
      Return
      End
C
C
      Subroutine Bl30( Bits, Qij, n20, MatDef )
C
C  Matrix made by matblas from blosum35.iij
C  Blocks Database = /data/blocks_5.0/blocks.dat
C
      Implicit None
C
      Integer       i, j, k, n20
      Integer       IQ(210)
      Real*4        Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits  /  '.142'  /
      Data  Def(1:40)  / 'Blosum 30, Entropy = 0.1424, Transition ' /,
     *      Def(41:60) / 'frequencies matrix. ' /
C
C
C ***   The blosum.qij matrix holds the 276 lower triangular elements of
C       the Blosum30 Qij, transition frequency, matrix.
C
C = A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
C
C  BLOSUM Clustered Target Frequencies=qij                                      
C  Blocks Database = /data/blocks_5.0/blocks.dat                                
C  Cluster Percentage: >= 30                                                    
C  amino acid order = ARNDCQEGHIL 
C ***   IQ values are q(i,j) * 10,000.  Below Diagonal
      Data  ( IQ(i), i = 1, 110 )  /
     *       96,   38,  109,   31,   19,   55,   43,   26,   27,   95,
     1       14,   11,   10,   10,   70,   31,   31,   14,   18,    7,
     2       39,   44,   31,   24,   37,   18,   28,   94,   52,   32,
     3       32,   35,   11,   20,   35,  173,   16,   14,   11,   11,
     4        4,   10,   18,   15,   60,   40,   22,   24,   18,   12,
     5       16,   23,   36,   12,   72,   56,   39,   32,   43,   23,
     6       27,   51,   51,   22,   66,  139,   44,   43,   27,   32,
     7       10,   21,   53,   39,   13,   26,   44,   63,   18,   12,
     8        9,    8,    4,    7,   12,   13,    8,   14,   27,   17,
     9       12,   27,   23,   17,   11,    8,   10,   16,   22,    8,
     A       27,   55,   23,    8,   77,   28,   21,   12,   21,    7 /
      Data  ( IQ(i), i = 111, 210 )  /
     *       15,   30,   30,   14,   17,   27,   29,    5,   11,   91,
     1       56,   35,   28,   34,   13,   21,   38,   51,   17,   33,
     2       47,   42,   10,   27,   24,   75,   40,   19,   24,   24,
     3        9,   16,   23,   28,   10,   27,   46,   25,   10,   16,
     4       20,   41,   46,    5,    7,    2,    4,    3,    4,    7,
     5       11,    2,    5,    9,    6,    2,    7,    4,    5,    3,
     6       27,   14,   22,    8,   15,    4,   11,   16,   16,   10,
     7       18,   45,   17,    7,   24,   11,   17,   14,    9,   44,
     8       56,   31,   22,   27,   12,   15,   26,   33,   12,   63,
     9       74,   30,   15,   32,   17,   36,   36,    5,   24,   83 /
C
      Bits = IBits
      MatDef = Def
      k = 0
      Do 150 i = 1, 20
         Do 140 j = 1, i
            k = k + 1
            Qij(i,j) = 0.0001 * IQ(k)
            If( j .lt. i )    Qij(j,i) = Qij(i,j)
  140       continue
  150    continue
      Return
      End
C
C
      Subroutine Bl40( Bits, Qij, n20, MatDef )
C
C  Matrix made by matblas from blosum40.qij
C  Blocks Database = /data/blocks_5.0/blocks.dat
C
      Implicit None
C
      Integer       i, j, k, n20
      Integer       IQ(210)
      Real*4        Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits  /  '.380'  /
      Data  Def(1:40)  / 'Blosum 40, Entropy = 0.3795, Transition ' /,
     *      Def(41:60) / 'frequencies matrix. ' /
C
C ***   The IQ matrix holds the 276 lower triangular elements of
C       the Blosum40.qij transition frequency matrix time 10,000
C
C = A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
C
C  BLOSUM Clustered Target Frequencies=qij                                      
C  Blocks Database = /data/blocks_5.0/blocks.dat                                
C  Cluster Percentage: >= 40                                                    
C  Amino Acid order = ARNDCQEGHILKMFPSTWYV
C ***   IQ values are q(i,j) * 10,000.  Below Diagonal
      Data  ( IQ(i), i = 1, 110 )  /
     *      148,   29,  109,   29,   19,   69,   32,   21,   31,  126,
     1       15,    7,    7,    9,   93,   26,   24,   17,   16,    4,
     2       48,   40,   26,   23,   48,   10,   32,  118,   66,   23,
     3       35,   31,   12,   20,   28,  260,   14,   12,   13,   14,
     4        3,   10,   15,   14,   60,   37,   17,   17,   16,    8,
     5       12,   18,   24,    9,  105,   50,   30,   21,   27,   15,
     6       23,   36,   34,   17,   82,  209,   41,   51,   27,   30,
     7       10,   26,   45,   35,   13,   23,   37,   99,   17,   10,
     8        7,    7,    3,    7,   10,   13,    7,   18,   37,   12,
     9       18,   22,   16,   13,   13,    8,    8,   15,   21,    9,
     A       31,   55,   18,   11,  105,   27,   14,   14,   19,    5 /
      Data  ( IQ(i), i = 111, 210 )  /
     *       13,   26,   28,    8,   20,   21,   23,    7,   10,  151,
     1       60,   25,   30,   30,   13,   26,   34,   52,   13,   25,
     2       34,   34,   11,   19,   22,   89,   40,   19,   22,   24,
     3       11,   15,   25,   29,    9,   28,   39,   29,   11,   19,
     4       22,   43,   70,    7,    5,    3,    3,    1,    4,    5,
     5        8,    2,    5,    9,    5,    2,    8,    3,    4,    3,
     6       45,   19,   16,   10,   11,    4,   10,   14,   17,   12,
     7       20,   31,   17,   10,   32,    9,   15,   14,    8,   60,
     8       54,   23,   17,   22,   12,   14,   25,   28,    8,   83,
     9       82,   27,   18,   31,   18,   32,   40,    5,   20,  113 /
C
      Bits = IBits
      MatDef = Def
      k = 0
      Do 150 i = 1, 20
         Do 140 j = 1, i
            k = k + 1
            Qij(i,j) = 0.0001 * IQ(k)
            If( j .lt. i )    Qij(j,i) = Qij(i,j)
  140       continue
  150    continue
      Return
      End
C
C
      Subroutine Bl50( Bits, Qij, n20, MatDef )
C
C  Matrix made by matblas from blosum50.qij
C  Blocks Database = /data/blocks_5.0/blocks.dat
C
      Implicit None
C
      Integer       i, j, k, n20
      Integer       IQ(210)
      Real*4        Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits  / '.481' /
      Data  Def(1:40)  / 'Blosum 50, Entropy = 0.4808, Transition ' /,
     *      Def(41:60) / 'frequencies matrix. ' /
C
C ***   The blosum matrix holds the 276 lower triangular elements of
C       the Blosum50 transitions frequencies matrix
C
C = A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
C
C  BLOSUM Clustered Target Frequencies=qij                                      
C  Blocks Database = /data/blocks_5.0/blocks.dat                                
C  Cluster Percentage: >= 50                                                    
C  Amino acid order = ARNDCQEGHILKMFPSTWYV      
C ***   IQ values are q(i,j) * 10,000.  Below Diagonal
      Data  ( IQ(i), i = 1, 110 )  /
     *      192,   27,  152,   24,   20,  101,   26,   19,   35,  161,
     1       15,    5,    6,    5,   91,   22,   25,   16,   17,    4,
     2       57,   34,   29,   23,   48,    6,   33,  141,   62,   20,
     3       31,   28,    9,   17,   23,  316,   12,   13,   15,   11,
     4        3,   10,   13,   11,   64,   35,   15,   13,   12,    8,
     5       11,   15,   18,    7,  140,   48,   28,   17,   18,   14,
     6       19,   26,   27,   13,  104,  304,   33,   64,   27,   26,
     7        6,   29,   43,   28,   14,   17,   27,  130,   16,    9,
     8        7,    5,    4,    8,    8,    9,    5,   22,   42,   10,
     9       29,   20,   12,    9,    8,    6,    7,   12,   15,    9,
     A       30,   58,   12,   12,  154,   22,   11,   11,   15,    4 /
      Data  ( IQ(i), i = 111, 210 )  /
     *       11,   19,   19,    6,   13,   17,   18,    5,    7,  171,
     1       62,   25,   32,   28,   11,   22,   30,   44,   12,   21,
     2       29,   31,   10,   16,   18,  111,   39,   21,   26,   22,
     3       10,   15,   24,   25,    9,   29,   38,   26,   11,   15,
     4       16,   47,  100,    5,    4,    2,    2,    1,    3,    4,
     5        5,    2,    5,    8,    4,    3,    9,    2,    3,    4,
     6       59,   15,   13,    9,    9,    4,    9,   12,   12,   13,
     7       18,   27,   13,    7,   39,    6,   13,   12,    8,   77,
     8       54,   20,   15,   16,   13,   14,   20,   22,    7,  107,
     9       92,   22,   21,   30,   16,   29,   41,    5,   18,  164 /
C
      Bits = IBits
      MatDef = Def
      k = 0
      Do 150 i = 1, 20
         Do 140 j = 1, i
            k = k + 1
            Qij(i,j) = 0.0001 * IQ(k)
            If( j .lt. i )    Qij(j,i) = Qij(i,j)
  140       continue
  150    continue
      Return
      End
C
C
      Subroutine Bl60( Bits, Qij, n20, MatDef )
C
C  Matrix made by matblas from blosum60.qij
C  Blocks Database = /data/blocks_5.0/blocks.dat
C
      Implicit None
C
      Integer       i, j, k, n20
      Integer       IQ(210)
      Real*4        Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits  / '.660' /
      Data  Def(1:40)  / 'Blosum 60, Entropy = 0.6603, Transition ' /,
     *      Def(41:60) / 'frequencies matrix. ' /
C
C ***   The blosum.qij matrix holds the 276 lower triangular elements of
C       the Blosum60 transitions frequencies matrix
C
C = A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
C
C  BLOSUM Clustered Target Frequencies=qij                                      
C  Blocks Database = /data/blocks_5.0/blocks.dat                                
C  Cluster Percentage: >= 60                                                    
C  Amino Acid order = ARNDCQEGHILKMFPSTWYV      
C ***   IQ values are q(i,j) * 10,000.  Below Diagonal
      Data  ( IQ(i), i = 1, 110 )  /
     *      210,   24,  174,   20,   19,  129,   22,   16,   37,  197,
     1       16,    4,    5,    4,  116,   20,   25,   15,   17,    3,
     2       69,   30,   27,   22,   50,    4,   35,  161,   58,   17,
     3       28,   25,    8,   13,   20,  364,   11,   13,   14,   10,
     4        2,   11,   13,   10,   84,   33,   13,   10,   13,   10,
     5        9,   13,   15,    6,  173,   46,   25,   14,   16,   15,
     6       17,   21,   21,   10,  111,  362,   33,   63,   24,   25,
     7        5,   31,   42,   25,   12,   16,   25,  152,   14,    8,
     8        5,    5,    4,    8,    7,    8,    4,   25,   48,    9,
     9       40,   17,   10,    8,    8,    5,    6,    9,   12,    8,
     A       31,   56,   10,   12,  182,   21,   10,    9,   13,    4 /
      Data  ( IQ(i), i = 111, 210 )  /
     *        9,   15,   14,    5,   10,   15,   17,    4,    6,  194,
     1       63,   23,   31,   28,   11,   19,   30,   40,   11,   18,
     2       25,   31,    9,   13,   17,  127,   38,   19,   23,   20,
     3       10,   14,   21,   22,    8,   27,   34,   24,   11,   12,
     4       14,   48,  123,    4,    3,    2,    2,    2,    3,    3,
     5        4,    2,    4,    8,    3,    2,    9,    2,    3,    3,
     6       66,   14,   10,    7,    6,    4,    7,    9,    9,   14,
     7       15,   23,   10,    6,   42,    5,   11,   10,    9,   99,
     8       52,   17,   12,   13,   14,   12,   17,   19,    6,  117,
     9       95,   20,   24,   26,   13,   25,   38,    4,   16,  192 /
C
      Bits = IBits
      MatDef = Def
      k = 0
      Do 150 i = 1, 20
         Do 140 j = 1, i
            k = k + 1
            Qij(i,j) = 0.0001 * IQ(k)
            If( j .lt. i )    Qij(j,i) = Qij(i,j)
  140       continue
  150    continue
      Return
      End
C
C
      Subroutine Bl70( Bits, Qij, n20, MatDef )
C
C  Matrix made by matblas from blosum70.qij
C  Blocks Database = /data/blocks_5.0/blocks.dat
C
      Implicit None
C
      Integer       i, j, k, n20
      Integer       IQ(210)
      Real*4        Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits  / '.839' /
      Data  Def(1:40)  / 'Blosum 70, Entropy = 0.8391, Transition ' /,
     *      Def(41:60) / 'frequencies matrix. ' /
C
C ***   The IQ matrix holds the 276 lower triangular elements of
C       the Blosum70 transitions frequencies matrix
C
C = A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
C
C  BLOSUM Clustered Target Frequencies=qij                                      
C  Blocks Database = /data/blocks_5.0/blocks.dat                                
C  Cluster Percentage: >= 70                                                    
C  Amino Acid Order = ARNDCQEGHILKMFPSTWYV      
C ***   IQ values are q(i,j) * 10,000.  Below Diagonal
      Data  ( IQ(i), i = 1, 110 )  /
     *      236,   22,  193,   18,   18,  155,   20,   14,   38,  238,
     1       16,    4,    4,    4,  144,   18,   25,   14,   15,    3,
     2       83,   28,   25,   21,   49,    3,   34,  183,   56,   16,
     3       27,   24,    7,   12,   18,  421,   10,   12,   13,    9,
     4        2,   11,   14,    8,   99,   29,   11,    9,   10,   11,
     5        8,   11,   11,    5,  202,   40,   20,   12,   12,   15,
     6       15,   17,   18,    9,  117,  417,   31,   62,   23,   22,
     7        4,   30,   39,   23,   11,   13,   23,  173,   12,    7,
     8        5,    4,    4,    7,    6,    6,    3,   25,   52,    8,
     9       46,   15,    8,    6,    7,    6,    5,    7,   11,    8,
     A       29,   55,    8,   11,  197,   22,    9,    7,   10,    3 /
      Data  ( IQ(i), i = 111, 210 )  /
     *        8,   13,   12,    4,    8,   13,   14,    4,    5,  206,
     1       64,   21,   29,   26,   10,   17,   28,   35,   11,   16,
     2       23,   28,    8,   11,   16,  142,   36,   16,   21,   18,
     3        9,   13,   20,   20,    7,   25,   30,   22,   10,   11,
     4       12,   46,  138,    4,    2,    2,    1,    2,    2,    2,
     5        4,    2,    4,    7,    2,    2,    9,    1,    3,    2,
     6       75,   12,    9,    6,    5,    3,    6,    7,    7,   17,
     7       14,   22,    9,    5,   43,    4,   10,    9,    9,  124,
     8       48,   15,   11,   11,   13,   11,   15,   16,    6,  119,
     9       94,   17,   22,   25,   11,   21,   34,    4,   14,  218 /
C
      Bits = IBits
      MatDef = Def
      k = 0
      Do 150 i = 1, 20
         Do 140 j = 1, i
            k = k + 1
            Qij(i,j) = 0.0001 * IQ(k)
            If( j .lt. i )    Qij(j,i) = Qij(i,j)
  140       continue
  150    continue
      Return
      End
C
C
      Subroutine Bl80( Bits, Qij, n20, MatDef )
C
C  Matrix made by matblas from blosum80.qij
C  Blocks Database = /data/blocks_5.0/blocks.dat
C
      Implicit None
C
      Integer       i, j, k, n20
      Integer       IQ(210)
      Real*4        Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits  / '.987' /
      Data  Def(1:40)  / 'Blosum 80, Entropy = 0.9868, Transition ' /,
     *      Def(41:60) / 'frequencies matrix. ' /
C
C ***   The IQ matrix holds the 276 lower triangular elements of
C       the Blosum70 transitions frequencies matrix
C
C = A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
C
C  BLOSUM Clustered Target Frequencies=qij                                      
C  Blocks Database = /data/blocks_5.0/blocks.dat                                
C  Cluster Percentage: >= 80                                                    
C  Amino Acid Order = ARNDCQEGHILKMFPSTWYV      
C ***   IQ values are q(i,j) * 10,000.  Below Diagonal
      Data  ( IQ(i), i = 1, 110 )  /
     *      252,   20,  210,   16,   17,  166,   18,   13,   37,  262,
     1       15,    3,    4,    3,  172,   17,   24,   14,   14,    3,
     2       94,   28,   23,   19,   48,    3,   35,  208,   53,   15,
     3       25,   22,    6,   11,   17,  463,    9,   12,   12,    8,
     4        2,   11,   12,    8,  104,   27,   10,    7,    8,   11,
     5        7,   10,    9,    4,  220,   36,   18,   11,   11,   14,
     6       14,   15,   16,    8,  111,  442,   29,   61,   22,   20,
     7        4,   28,   36,   20,   10,   12,   19,  190,   11,    6,
     8        4,    3,    4,    7,    6,    5,    3,   25,   52,    7,
     9       53,   14,    7,    6,    6,    5,    5,    6,    9,    7,
     A       27,   52,    7,   10,  211,   21,    9,    7,    9,    3 /
      Data  ( IQ(i), i = 111, 210 )  /
     *        7,   12,   10,    4,    7,   12,   12,    3,    4,  221,
     1       64,   20,   29,   24,   10,   17,   26,   34,   10,   15,
     2       21,   26,    7,   10,   14,  167,   36,   15,   20,   16,
     3        9,   12,   19,   19,    7,   24,   28,   20,    9,   11,
     4       11,   48,  156,    3,    2,    1,    1,    1,    2,    2,
     5        3,    1,    3,    6,    2,    2,    7,    1,    2,    2,
     6       87,   11,    7,    6,    5,    3,    5,    6,    6,   16,
     7       13,   20,    8,    5,   46,    3,    9,    8,   10,  148,
     8       46,   13,    9,   10,   13,   10,   15,   14,    5,  123,
     9       89,   15,   22,   22,   10,   21,   33,    4,   12,  246 /
C    
      Bits = IBits
      MatDef = Def
      k = 0
      Do 150 i = 1, 20
         Do 140 j = 1, i
            k = k + 1
            Qij(i,j) = 0.0001 * IQ(k)
            If( j .lt. i )    Qij(j,i) = Qij(i,j)
  140       continue
  150    continue
      Return
      End
C
C
      Subroutine Bl90( Bits, Qij, n20, MatDef )
C
C  Matrix made by matblas from blosum90.qij
C  Blocks Database = /data/blocks_5.0/blocks.dat
C
      Implicit None
C
      Integer       i, j, k, n20
      Integer       IQ(210)
      Real*4        Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits  / '1.18' /
      Data  Def(1:40)  / 'Blosum 90, Entropy = 1.1806, Transition ' /,
     *      Def(41:60) / 'frequencies matrix. ' /
C
C ***   The blosum matrix holds the 276 lower triangular elements of
C       the Blosum90.qij transition frequencies matrix
C
C = A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
C
C  BLOSUM Clustered Target Frequencies=qij                                      
C  Blocks Database = /data/blocks_5.0/blocks.dat                                
C  Cluster Percentage: >= 90                                                    
C  Amino Acid Order = ARNDCQEGHILKMFPSTWYV      
C ***   IQ values are q(i,j) * 10,000.  Below Diagonal
      Data  ( IQ(i), i = 1, 110 )  /
     *      289,   19,  229,   14,   15,  178,   16,   11,   34,  289,
     1       14,    3,    4,    3,  190,   16,   21,   13,   12,    2,
     2      113,   27,   20,   17,   48,    2,   34,  241,   49,   13,
     3       22,   21,    5,   10,   15,  521,    8,   11,   11,    7,
     4        1,   10,   11,    7,  115,   23,    9,    6,    6,   10,
     5        6,    9,    8,    4,  245,   31,   16,    9,    9,   12,
     6       12,   13,   14,    7,  102,  473,   26,   57,   20,   18,
     7        3,   27,   32,   18,    9,   10,   16,  211,   10,    6,
     8        3,    3,    3,    7,    5,    5,    2,   24,   49,    7,
     9       68,   12,    6,    5,    5,    5,    4,    5,    7,    6,
     A       24,   47,    7,    9,  228,   20,    8,    6,    8,    3 /
      Data  ( IQ(i), i = 111, 210 )  /
     *        6,   11,    9,    4,    7,   10,   11,    3,    4,  240,
     1       61,   18,   26,   22,    9,   15,   24,   32,    9,   13,
     2       19,   23,    6,    9,   13,  202,   35,   13,   18,   15,
     3        8,   12,   17,   16,    6,   21,   25,   19,    9,    9,
     4       10,   49,  193,    3,    2,    1,    1,    1,    2,    2,
     5        3,    1,    3,    6,    1,    2,    7,    1,    2,    2,
     6      101,   10,    6,    5,    4,    3,    5,    5,    5,   14,
     7       11,   18,    7,    4,   46,    3,    9,    8,   10,  173,
     8       42,   11,    8,    8,   13,    9,   14,   12,    5,  120,
     9       81,   13,   20,   20,    9,   19,   30,    3,   11,  287 /
C
      Bits = IBits
      MatDef = Def
      k = 0
      Do 150 i = 1, 20
         Do 140 j = 1, i
            k = k + 1
            Qij(i,j) = 0.0001 * IQ(k)
            If( j .lt. i )    Qij(j,i) = Qij(i,j)
  140       continue
  150    continue
      Return
      End
C
C
      Subroutine Bl100( Bits, Qij, n20, MatDef )
C
C  Matrix made by matblas from blosum100_3.iij
C  Blocks Database = /data/blocks_5.0/blocks.dat
C
      Implicit None
C
      Integer       i, j, k, n20
      Integer       IQ(210)
      Real*4        Qij(n20,n20)
      Character*4   IBits, Bits
      Character*60  MatDef, Def
C
      Data  IBits  / '1.45' /
      Data  Def(1:40)  / 'Blosum 90, Entropy = 1.4516, Transition ' /,
     *      Def(41:60) / 'frequencies matrix. ' /
C
C ***   The blosum matrix holds the 276 lower triangular elements of
C       the Blosum100.qij transition frequencies matrix
C
C = A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
C
C  BLOSUM Clustered Target Frequencies=qij                                      
C  Blocks Database = /data/blocks_5.0/blocks.dat                                
C  Cluster Percentage: >= 100                                                   
C  A      R      N      D      C      Q      E      G      H      I      L      
C ***   IQ values are q(i,j) * 10,000.  Below Diagonal
      Data  ( IQ(i), i = 1, 110 )  /
     *      340,   16,  266,   13,   13,  204,   14,    9,   30,  321,
     1       13,    2,    4,    2,  201,   15,   17,   11,   11,    2,
     2      139,   22,   17,   15,   47,    2,   29,  297,   42,   11,
     3       19,   18,    4,    8,   12,  581,    7,   10,   10,    6,
     4        1,   10,    9,    5,  138,   19,    7,    5,    5,    9,
     5        5,    7,    6,    3,  284,   27,   13,    7,    7,   10,
     6       11,   10,   12,    6,   90,  526,   22,   52,   17,   15,
     7        3,   24,   27,   15,    7,    8,   13,  240,    8,    5,
     8        3,    2,    3,    6,    4,    4,    2,   22,   45,    6,
     9       93,   10,    5,    4,    4,    5,    4,    4,    6,    5,
     A       20,   41,    6,    9,  257,   17,    6,    4,    7,    2 /
      Data  ( IQ(i), i = 111, 210 )  /
     *        6,    9,    8,    3,    5,    8,    9,    3,    3,  263,
     1       58,   15,   23,   19,    9,   13,   22,   29,    8,   11,
     2       15,   19,    5,    8,   11,  240,   32,   12,   17,   13,
     3        8,   10,   15,   15,    6,   19,   21,   15,    8,    8,
     4        8,   47,  243,    3,    2,    1,    1,    1,    2,    1,
     5        2,    1,    2,    4,    1,    2,    7,    1,    2,    2,
     6      107,    8,    5,    5,    4,    3,    4,    4,    4,   12,
     7       10,   14,    5,    3,   44,    2,    7,    6,    9,  197,
     8       37,    9,    7,    7,   11,    7,   12,    9,    4,  112,
     9       71,   11,   19,   17,    7,   17,   29,    3,    9,  329 /
C
      Bits = IBits
      MatDef = Def
      k = 0
      Do 150 i = 1, 20
         Do 140 j = 1, i
            k = k + 1
            Qij(i,j) = 0.0001 * IQ(k)
            If( j .lt. i )    Qij(j,i) = Qij(i,j)
  140       continue
  150    continue
      Return
      End
C
C
      Subroutine pam80( Bits, Qij, n20, MatDef )
C
      INTEGER       i, j, k, n20
      Integer       IQ( 400 )
      Real*4        Qij( n20, n20 )
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits   /  '1.68'  /
      Data  Def(1:41)  / 'Dayhoff PAM 80 Transiton Frequency Matrix' /,
     *      Def(42:60) / ':  Entropy = 1.6815' /
C
C   Pam 80 Amino Acid Transition Frequency matrix
C
C ***   Order of amino acdids in PAM matrices
C      Character*1  AA(20)
C      DATA  AA  / 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
C     *            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' /
C
C       PAM Transition Frequency Matrix, Dist=80 Entropy=1.6815
C ***   IQ values are q(i,j) * 10,000.  Full Matrix
C
      Data  ( IQ(i), i = 1, 100 )  /
     *     3848,  218,  571,  579,  231,  454,  714,  954,  222,  436,
     1      265,  246,  353,  153,  976, 1215, 1211,   46,  144,  751,
     2      100, 5143,  170,   69,   72,  446,   89,   41,  471,  152,
     3       81,  850,  235,   68,  227,  255,  131,  431,   28,   84,
     4      267,  172, 2719,  979,   48,  280,  426,  303,  685,  141,
     5       70,  483,   86,   78,  173,  572,  385,   66,  173,  103,
     6      332,   78, 1138, 3753,   25,  467, 1595,  339,  307,  102,
     7       39,  242,   56,   24,  140,  328,  236,   17,   45,  103,
     8       82,   68,   38,   16, 8072,   17,   16,   24,   65,   75,
     9       13,   17,   19,   22,   80,  224,   90,    9,  181,  118 /
      Data  ( IQ(i), i = 101, 200 )  /
     *      184,  405,  269,  383,   19, 3926,  884,   99,  872,   92,
     1      154,  306,  190,   28,  289,  155,  141,   26,   32,   89,
     2      419,  106,  515, 1684,   25, 1148, 3884,  274,  265,  158,
     3       87,  252,  106,   27,  219,  271,  204,   13,   63,  147,
     4      961,  138,  629,  619,  123,  261,  481, 6122,  154,  123,
     5      105,  199,  130,   87,  318,  927,  393,   30,   45,  321,
     6       91,  378,  586,  224,   64,  757,  168,   42, 5082,   43,
     7       74,  136,   45,  122,  177,  119,   96,   75,  212,   66,
     8      165,  113,  133,   83,  117,   87,  106,   46,   52, 3886,
     9      457,  123,  527,  326,   52,  108,  299,   19,  101, 1159 /
      Data  ( IQ(i), i = 201, 300 )  /
     *      220,  122,  181,   68,   33,  320,  114,   97,  254, 1110,
     1     6711,  178, 1940,  752,  200,  132,  242,  256,  209,  854,
     2      232, 1654,  937,  444,   40,  625,  428,  187,  304,  254,
     3      124, 5770,  826,   46,  252,  471,  549,   78,   86,  142,
     4       59,   70,   34,   17,    9,   86,   23,   15,   23,  220,
     5      344,  165, 3719,   69,   22,   54,   95,   11,   15,  187,
     6       75,   65,   81,   22,   30,   29,   23,   67,  131,  369,
     7      349,   24,  235, 6641,   22,  106,   85,  201, 1477,   81,
     8      573,  286,  200,  144,   95,  389,  220,  192,  287,  112,
     9      135,  164,  112,   77, 5637,  522,  305,   26,   31,  161 /
      Data  ( IQ(i), i = 301, 400 )  /
     *      973,  464,  969,  471,  490,  298,  385,  714,  243,  211,
     1      115,  413,  235,  162,  736, 3196, 1217,  237,  145,  240,
     2      827,  207,  554,  295,  129,  217,  231,  261,  143,  473,
     3      162,  402,  315,  102,  368, 1026, 3859,   40,  128,  451,
     4        6,  108,    7,    2,    4,    5,    2,    4,    6,    4,
     5        3,   10,    4,   65,    5,   45,    8, 8257,   66,    2,
     6       64,   19,  132,   33,  181,   29,   55,   19,  209,   93,
     7       89,   16,   38, 1109,   17,   71,   69,  137, 6566,   67,
     8      538,  145,  150,  121,  193,  163,  162,  204,  170, 2007,
     9      628,  125,  792,  166,  214,  222,  498,   23,  135, 4885 /
C
C
C ***   Code to put PAM80 qij values into square scoring matrix, Qij.
C
      MatDef = Def
      Bits = IBits
      k = 1
      Do 110 j = 1, 20
         Do 100 i = 1, 20, 1
            Qij(i,j) = 0.0001 * Float( IQ( k ) ) 
            k = k + 1
  100       continue
  110    continue
      RETURN
      END
C
C
      Subroutine pam120( Bits, Qij, n20, MatDef )
C
      INTEGER       i, j, k, n20
      Integer       IQ( 400 )
      Real*4        Qij( n20, n20 )
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits   /  '1.20'  /
      Data  Def(1:40)  / 'Dayhoff PAM 120 Transiton Frequency Matr' /,
     *      Def(41:60) / 'ix: Entropy = 1.2017' /
C
C   Pam 120 Amino Acid Transition Frequency matrix
C
C ***   Order of amino acdids in PAM matrices
C      Character*1  AA(20)
C      DATA  AA  / 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
C     *            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' /
C
C       PAM Transition Frequency Matrix, Dist=120 Entropy=1.2017
C ***   IQ values are q(i,j) * 10,000.  Full Matrix.
      Data  ( IQ(i), i = 1, 100 )  /
     *     2658,  331,  727,  735,  327,  583,  830, 1113,  343,  580,
     1      365,  378,  463,  223, 1124, 1279, 1302,   86,  206,  859,
     2      151, 3796,  254,  130,  101,  517,  154,   78,  553,  193,
     3      119,  968,  300,   96,  291,  302,  197,  538,   54,  125,
     4      343,  257, 1626,  928,   83,  371,  527,  372,  700,  180,
     5      104,  529,  143,  112,  244,  568,  438,   92,  205,  155,
     6      417,  147, 1078, 2578,   51,  598, 1547,  432,  409,  153,
     7       74,  331,  105,   47,  221,  416,  323,   33,   78,  158,
     8      118,   95,   65,   33, 7262,   34,   32,   45,   90,  106,
     9       26,   34,   36,   43,  116,  266,  130,   18,  238,  155 /
      Data  ( IQ(i), i = 101, 200 )  /
     *      242,  471,  353,  489,   37, 2605,  892,  150,  926,  135,
     1      193,  377,  232,   54,  349,  216,  198,   49,   59,  132,
     2      486,  185,  640, 1632,   50, 1161, 2691,  369,  386,  208,
     3      130,  337,  162,   53,  304,  363,  296,   28,   91,  206,
     4     1124,  227,  779,  788,  197,  387,  648, 4904,  255,  223,
     5      168,  308,  212,  134,  475, 1076,  580,   60,   88,  434,
     6      135,  446,  597,  296,   88,  804,  251,   79, 3705,   77,
     7      104,  204,   83,  162,  230,  173,  143,  108,  263,   95,
     8      225,  149,  166,  119,  155,  128,  140,   85,   92, 2611,
     9      555,  165,  603,  389,   96,  158,  345,   38,  151, 1186 /
      Data  ( IQ(i), i = 201, 300 )  /
     *      313,  197,  253,  129,   68,  407,  181,  155,  343, 1344,
     1     5601,  265, 2181,  966,  278,  212,  348,  349,  327, 1078,
     2      352, 1881, 1030,  596,   80,  774,  567,  284,  465,  343,
     3      198, 4531,  923,   89,  370,  601,  671,  145,  131,  229,
     4       77,   95,   55,   32,   17,  100,   41,   29,   41,  252,
     5      386,  184, 2331,   94,   39,   71,  114,   21,   29,  220,
     6      109,   91,  116,   43,   58,   55,   45,   94,  178,  439,
     7      450,   47,  307, 5509,   44,  137,  124,  277, 1806,  141,
     8      658,  367,  291,  229,  143,  471,  305,  283,  369,  175,
     9      188,  240,  174,  113, 4307,  603,  411,   51,   60,  231 /
      Data  ( IQ(i), i = 301, 400 )  /
     *     1024,  546,  962,  602,  580,  412,  511,  832,  354,  309,
     1      181,  527,  314,  213,  847, 2086, 1225,  291,  203,  351,
     2      886,  303,  629,  401,  197,  306,  339,  384,  225,  548,
     3      234,  492,  394,  156,  491, 1033, 2589,   72,  175,  546,
     4       11,  134,   13,    5,    7,   10,    5,    7,   12,    8,
     5        6,   18,    8,   89,   10,   54,   14, 7507,   90,    4,
     6       88,   37,  157,   56,  239,   52,   74,   37,  256,  138,
     7      134,   31,   70, 1358,   34,   98,   96,  190, 5404,   96,
     8      612,  207,  233,  192,  262,  236,  234,  281,  229, 2058,
     9      794,  198,  927,  264,  298,  321,  600,   47,  193, 3620 /
C
C
C ***   Code to put PAM120 qij values into square scoring matrix, Qij.
C
      MatDef = Def
      Bits = IBits
      k = 1
      Do 110 j = 1, 20
         Do 100 i = 1, 20, 1
            Qij(i,j) = 0.0001 * Float( IQ( k ) ) 
            k = k + 1
  100       continue
  110    continue
      RETURN
      END
C
C
      Subroutine pam160( Bits, Qij, n20, MatDef )
C
      INTEGER       i, j, k, n20
      Integer       IQ( 400 )
      Real*4        Qij( n20, n20 )
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits   /  '.899'  /
      Data  Def(1:40)  / 'Dayhoff PAM 120 Transiton Frequency Matr' /,
     *      Def(41:60) / 'ix:  Entropy = .8991' /
C
C   Pam 80 Amino Acid Transition Frequency matrix
C
C ***   Order of amino acdids in PAM matrices
C      Character*1  AA(20)
C      DATA  AA  / 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
C     *            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' /
C
C       PAM Qij Transition Frequency matrix, Dist=160 Entropy=0.8991
C ***   IQ values are q(i,j) * 10,000.  Full Matrix.
      Data  ( IQ(i), i = 1, 100 )  /
     *     1989,  433,  821,  834,  406,  674,  889, 1178,  450,  679,
     1      449,  492,  546,  289, 1178, 1252, 1288,  130,  264,  904,
     2      197, 2861,  318,  190,  127,  545,  214,  119,  587,  223,
     3      153,  992,  341,  122,  334,  332,  251,  603,   82,  163,
     4      389,  321, 1091,  833,  116,  429,  570,  413,  667,  213,
     5      136,  538,  191,  142,  301,  539,  458,  116,  224,  200,
     6      469,  217,  966, 1892,   80,  664, 1389,  491,  476,  198,
     7      112,  396,  156,   75,  293,  469,  388,   52,  109,  208,
     8      148,  118,   90,   52, 6539,   53,   50,   68,  112,  134,
     9       42,   54,   55,   66,  147,  290,  163,   29,  280,  184 /
      Data  ( IQ(i), i = 101, 200 )  /
     *      283,  498,  407,  541,   58, 1811,  833,  198,  895,  174,
     1      220,  420,  260,   81,  382,  265,  245,   73,   87,  171,
     2      520,  258,  694, 1466,   80, 1085, 1986,  438,  474,  250,
     3      171,  403,  212,   83,  371,  429,  369,   48,  119,  256,
     4     1193,  317,  871,  899,  269,  498,  771, 3992,  353,  323,
     5      235,  409,  293,  183,  609, 1140,  723,   94,  137,  525,
     6      175,  474,  566,  343,  110,  776,  311,  116, 2744,  111,
     7      132,  260,  120,  192,  268,  214,  183,  136,  292,  122,
     8      267,  178,  194,  151,  184,  165,  169,  124,  130, 1854,
     9      606,  198,  626,  423,  138,  199,  365,   59,  196, 1108 /
      Data  ( IQ(i), i = 201, 300 )  /
     *      395,  272,  319,  195,  109,  471,  249,  216,  415, 1463,
     1     4730,  344, 2218, 1111,  348,  290,  440,  425,  439, 1213,
     2      456, 1927, 1049,  703,  126,  866,  670,  374,  597,  416,
     3      270, 3636,  943,  139,  473,  687,  744,  215,  176,  314,
     4       91,  112,   73,   48,   25,  109,   57,   43,   59,  262,
     5      392,  188, 1502,  115,   57,   84,  126,   31,   44,  232,
     6      140,  115,  148,   68,   89,   84,   70,  119,  217,  477,
     7      518,   73,  358, 4624,   69,  164,  158,  339, 1973,  198,
     8      689,  422,  363,  306,  188,  517,  373,  359,  426,  234,
     9      234,  306,  231,  148, 3335,  637,  484,   78,   93,  290 /
      Data  ( IQ(i), i = 301, 400 )  /
     *     1003,  595,  913,  680,  631,  501,  600,  883,  441,  390,
     1      245,  602,  377,  257,  892, 1515, 1148,  328,  252,  438,
     2      874,  380,  655,  478,  256,  379,  425,  478,  299,  583,
     3      297,  547,  447,  205,  573,  969, 1844,  106,  215,  596,
     4       16,  150,   19,    9,   11,   16,    9,   11,   19,   13,
     5       10,   26,   13,  107,   15,   59,   20, 6826,  109,    8,
     6      109,   58,  173,   78,  280,   76,   91,   57,  283,  177,
     7      175,   49,  104, 1484,   53,  121,  120,  234, 4493,  122,
     8      643,  264,  306,  258,  318,  300,  296,  345,  281, 1926,
     9      894,  267,  980,  351,  368,  398,  653,   77,  246, 2781 /
C
C
C ***   Code to put PAM160 qij values into square scoring matrix, Qij.
C
      MatDef = Def
      Bits = IBits
      k = 1
      Do 110 j = 1, 20
         Do 100 i = 1, 20, 1
            Qij(i,j) = 0.0001 * Float( IQ( k ) ) 
            k = k + 1
  100       continue
  110    continue
      RETURN
      END
C
C
      Subroutine pam200( Bits, Qij, n20, MatDef )
C
      INTEGER       i, j, k, n20
      Integer       IQ( 400 )
      Real*4        Qij( n20, n20 )
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits   /  '.697'  /
      Data  Def(1:40)  / 'Dayhoff PAM 200 Transiton Frequency Matr' /,
     *      Def(41:60) / 'ix:  Entropy = .6969' /
C
C   Pam 200 Amino Acid Transition Frequency matrix
C
C ***   Order of amino acdids in PAM matrices
C      Character*1  AA(20)
C      DATA  AA  / 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
C     *            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' /
C
C       PAM Transition Frequency Matrix, Dist=200 Entropy=0.6969
C ***   IQ values are q(i,j) * 10,000.  Full Matrix
C
      Data  ( IQ(i), i = 1, 100 )  /
     *     1603,  520,  876,  894,  473,  740,  922, 1194,  542,  744,
     1      519,  587,  611,  349, 1184, 1199, 1236,  174,  317,  919,
     2      237, 2207,  363,  244,  150,  550,  265,  158,  593,  248,
     3      185,  966,  366,  146,  364,  352,  294,  640,  111,  196,
     4      415,  366,  818,  743,  145,  463,  579,  437,  624,  242,
     5      168,  533,  230,  169,  344,  511,  462,  138,  238,  238,
     6      500,  280,  861, 1469,  111,  687, 1218,  529,  516,  238,
     7      150,  441,  203,  104,  353,  499,  433,   74,  137,  252,
     8      173,  137,  114,   73, 5892,   73,   70,   91,  132,  158,
     9       60,   74,   76,   91,  174,  304,  190,   41,  310,  207 /
      Data  ( IQ(i), i = 101, 200 )  /
     *      314,  505,  437,  558,   81, 1324,  758,  239,  831,  207,
     1      241,  444,  280,  109,  401,  303,  283,   98,  113,  206,
     2      538,  321,  706, 1285,  111,  990, 1546,  487,  531,  286,
     3      209,  451,  257,  116,  423,  474,  423,   71,  146,  298,
     4     1212,  402,  927,  970,  336,  593,  859, 3304,  445,  414,
     5      301,  502,  371,  232,  717, 1159,  825,  131,  187,  599,
     6      208,  480,  526,  370,  130,  721,  350,  151, 2067,  142,
     7      156,  301,  153,  216,  294,  245,  216,  161,  308,  148,
     8      295,  203,  218,  179,  208,  196,  195,  160,  165, 1390,
     9      626,  225,  622,  442,  177,  232,  373,   81,  234,  998 /
      Data  ( IQ(i), i = 201, 300 )  /
     *      465,  342,  381,  262,  155,  520,  314,  276,  474, 1510,
     1     4039,  414, 2151, 1204,  410,  362,  517,  488,  541, 1284,
     2      541, 1877, 1040,  775,  175,  920,  745,  454,  698,  477,
     3      336, 2981,  930,  191,  560,  743,  789,  281,  220,  390,
     4      102,  123,   88,   63,   34,  115,   71,   56,   74,  260,
     5      379,  185, 1004,  130,   72,   94,  132,   42,   59,  235,
     6      168,  138,  176,   94,  121,  113,   97,  141,  249,  498,
     7      563,   99,  395, 3924,   96,  187,  189,  389, 2031,  247,
     8      692,  461,  419,  370,  228,  544,  426,  419,  465,  286,
     9      273,  362,  280,  180, 2622,  647,  531,  107,  127,  339 /
      Data  ( IQ(i), i = 301, 400 )  /
     *      960,  627,  865,  725,  660,  569,  660,  900,  507,  453,
     1      303,  652,  430,  296,  902, 1209, 1059,  357,  294,  503,
     2      837,  439,  660,  532,  306,  439,  489,  545,  363,  598,
     3      349,  580,  482,  249,  625,  892, 1398,  140,  250,  620,
     4       21,  158,   24,   13,   15,   21,   13,   15,   26,   17,
     5       14,   34,   18,  121,   21,   62,   25, 6209,  125,   12,
     6      128,   78,  185,   98,  311,   98,  107,   76,  297,  210,
     7      211,   67,  137, 1529,   73,  141,  140,  271, 3774,  147,
     8      653,  315,  367,  318,  364,  355,  350,  398,  328, 1737,
     9      948,  329,  989,  425,  424,  455,  678,  110,  296, 2209 /
C
C
C ***   Code to put PAM200 qij values into square scoring matrix, Qij.
C
      MatDef = Def
      Bits = IBits
      k = 1
      Do 110 j = 1, 20
         Do 100 i = 1, 20, 1
            Qij(i,j) = 0.0001 * Float( IQ( k ) ) 
            k = k + 1
  100       continue
  110    continue
      RETURN
      END
C
C
      Subroutine pam250( Bits, Qij, n20, MatDef )
C
      INTEGER       i, j, k, n20
      Integer       IQ( 400 )
      Real*4        Qij( n20, n20 )
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits   /  '.528'  /
      Data  Def(1:40)  / 'Dayhoff PAM 250 Transiton Frequency Matr' /,
     *      Def(41:60) / 'ix:  Entropy = .5277' /
C
C   Pam 80 Amino Acid Transition Frequency matrix
C
C ***   Order of amino acdids in PAM matrices
C      Character*1  AA(20)
C      DATA  AA  / 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
C     *            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' /
C
C       PAM Transition Frequency Matrix, Dist=250 Entropy=0.5277
C ***   IQ values are q(i,j) * 10,000.  Full Matrix.
C
      Data  ( IQ(i), i = 1, 100 )  /
     *     1332,  609,  911,  935,  541,  797,  942, 1179,  633,  794,
     1      590,  680,  673,  417, 1162, 1131, 1163,  230,  378,  921,
     2      278, 1649,  399,  299,  177,  542,  316,  205,  581,  274,
     3      219,  902,  384,  175,  388,  370,  334,  660,  144,  233,
     4      431,  403,  648,  656,  179,  483,  567,  452,  573,  273,
     5      204,  520,  269,  199,  382,  485,  458,  163,  251,  276,
     6      520,  344,  760, 1142,  149,  682, 1033,  554,  541,  282,
     7      195,  479,  254,  141,  411,  518,  469,  102,  170,  299,
     8      198,  159,  140,   98, 5179,   99,   96,  119,  154,  183,
     9       83,  100,  100,  120,  202,  315,  216,   56,  336,  229 /
      Data  ( IQ(i), i = 101, 200 )  /
     *      341,  499,  453,  552,  109,  963,  671,  281,  741,  242,
     1      261,  459,  299,  141,  413,  336,  319,  127,  143,  242,
     2      549,  384,  693, 1090,  150,  879, 1204,  528,  570,  325,
     3      251,  494,  304,  156,  471,  509,  469,  102,  179,  341,
     4     1199,  500,  966, 1020,  413,  689,  930, 2670,  545,  513,
     5      379,  601,  459,  293,  819, 1151,  907,  179,  251,  670,
     6      241,  471,  480,  386,  152,  644,  377,  190, 1490,  176,
     7      182,  336,  189,  239,  316,  272,  249,  187,  317,  178,
     8      316,  229,  243,  210,  232,  228,  222,  199,  202, 1040,
     9      627,  252,  599,  454,  219,  262,  376,  108,  272,  861 /
      Data  ( IQ(i), i = 201, 300 )  /
     *      538,  421,  451,  341,  214,  569,  388,  349,  534, 1507,
     1     3366,  490, 2005, 1271,  479,  442,  594,  554,  650, 1315,
     2      624, 1753, 1015,  832,  236,  954,  808,  541,  785,  539,
     3      410, 2396,  898,  255,  648,  787,  820,  357,  273,  472,
     4      112,  133,  102,   79,   46,  120,   86,   72,   91,  250,
     5      353,  178,  647,  144,   89,  104,  138,   55,   76,  229,
     6      200,  165,  206,  126,  160,  148,  130,  167,  281,  510,
     7      596,  133,  426, 3244,  129,  213,  222,  438, 2006,  298,
     8      679,  492,  469,  433,  272,  560,  474,  474,  497,  341,
     9      315,  419,  333,  218, 1988,  642,  564,  142,  168,  388 /
      Data  ( IQ(i), i = 301, 400 )  /
     *      906,  655,  822,  754,  682,  628,  707,  895,  568,  513,
     1      368,  690,  483,  342,  892, 1009,  962,  387,  339,  561,
     2      785,  493,  655,  575,  356,  495,  544,  598,  426,  604,
     3      403,  604,  513,  298,  659,  810, 1074,  180,  289,  630,
     4       27,  162,   31,   19,   20,   28,   19,   21,   33,   24,
     5       20,   43,   24,  135,   27,   66,   32, 5518,  139,   17,
     6      149,  103,  196,  120,  337,  123,  127,  101,  304,  243,
     7      250,   91,  174, 1511,   97,  162,  163,  307, 3075,  175,
     8      654,  370,  427,  381,  410,  413,  406,  450,  379, 1502,
     9      972,  394,  967,  498,  479,  506,  687,  154,  352, 1729 /
C
C
C ***   Code to put PAM250 qij values into square scoring matrix, Qij.
C
      MatDef = Def
      Bits = IBits
      k = 1
      Do 110 j = 1, 20
         Do 100 i = 1, 20, 1
            Qij(i,j) = 0.0001 * Float( IQ( k ) ) 
            k = k + 1
  100       continue
  110    continue
      RETURN
      END
C
C
      Subroutine pam300( Bits, Qij, n20, MatDef )
C
      INTEGER       i, j, k, n20
      Integer       IQ( 400 )
      Real*4        Qij( n20, n20 )
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  IBits   /  '.415'  /
      Data  Def(1:40)  / 'Dayhoff PAM 300 Transiton Frequency Matr' /,
     *      Def(41:60) / 'ix:  Entropy = .4153' /
C
C   Pam 80 Amino Acid Transition Frequency matrix
C
C ***   Order of amino acdids in PAM matrices
C      Character*1  AA(20)
C      DATA  AA  / 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
C     *            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' /
C
C       PAM Transition Frequency Matrix, Dist=300 Entropy=0.4153
C ***   IQ values are q(i,j) * 10,000.  Full Matrix.
C
      Data  ( IQ(i), i = 1, 100 )  /
     *     1183,  680,  926,  953,  595,  836,  951, 1148,  703,  824,
     1      647,  750,  720,  477, 1127, 1076, 1101,  283,  433,  916,
     2      311, 1280,  420,  341,  200,  527,  353,  246,  560,  296,
     3      248,  828,  392,  201,  404,  383,  362,  662,  175,  265,
     4      438,  424,  564,  594,  208,  487,  546,  458,  533,  299,
     5      236,  506,  298,  225,  407,  468,  453,  186,  263,  306,
     6      528,  394,  688,  941,  186,  660,  891,  565,  550,  318,
     7      236,  502,  297,  177,  452,  525,  490,  131,  199,  336,
     8      219,  177,  163,  123, 4558,  123,  120,  144,  174,  204,
     9      106,  124,  124,  148,  224,  321,  236,   72,  353,  246 /
      Data  ( IQ(i), i = 101, 200 )  /
     *      359,  487,  456,  532,  137,  756,  602,  314,  660,  270,
     1      277,  463,  314,  171,  418,  359,  346,  153,  170,  271,
     2      553,  431,  669,  940,  189,  790,  992,  551,  585,  358,
     3      289,  521,  343,  195,  503,  528,  499,  133,  210,  375,
     4     1170,  584,  984, 1043,  480,  762,  972, 2216,  628,  593,
     5      452,  683,  535,  352,  889, 1128,  954,  229,  312,  725,
     6      266,  455,  443,  390,  172,  574,  387,  223, 1111,  205,
     7      204,  357,  219,  256,  329,  292,  273,  208,  319,  204,
     8      330,  251,  265,  236,  251,  255,  246,  231,  233,  830,
     9      612,  274,  569,  457,  252,  285,  374,  135,  302,  745 /
      Data  ( IQ(i), i = 201, 300 )  /
     *      598,  490,  513,  414,  274,  609,  455,  418,  584, 1468,
     1     2850,  553, 1842, 1300,  538,  511,  654,  608,  740, 1306,
     2      686, 1608,  989,  865,  295,  966,  847,  611,  839,  589,
     3      473, 1986,  867,  317,  715,  812,  838,  423,  323,  541,
     4      120,  138,  112,   93,   56,  125,   99,   85,  104,  238,
     5      324,  170,  452,  154,  102,  113,  141,   66,   91,  221,
     6      228,  190,  233,  158,  196,  181,  162,  192,  305,  515,
     7      610,  165,  446, 2723,  162,  236,  251,  475, 1917,  336,
     8      659,  513,  503,  478,  310,  568,  508,  513,  516,  385,
     9      351,  462,  375,  254, 1554,  631,  579,  176,  207,  424 /
      Data  ( IQ(i), i = 301, 400 )  /
     *      861,  674,  794,  766,  694,  668,  733,  878,  611,  557,
     1      424,  713,  525,  383,  873,  905,  891,  413,  379,  600,
     2      741,  530,  647,  601,  396,  536,  579,  628,  474,  604,
     3      445,  618,  533,  341,  673,  749,  893,  218,  323,  630,
     4       33,  162,   37,   25,   25,   34,   24,   26,   41,   30,
     5       26,   51,   30,  144,   33,   68,   38, 4905,  149,   22,
     6      168,  126,  206,  140,  353,  145,  145,  123,  306,  268,
     7      281,  114,  205, 1446,  120,  180,  182,  336, 2541,  199,
     8      650,  418,  474,  433,  446,  459,  452,  491,  424, 1304,
     9      966,  447,  931,  552,  521,  540,  685,  197,  402, 1413 /
C
C
C ***   Code to put PAM300 qij values into square scoring matrix, Qij.
C
      MatDef = Def
      Bits = IBits
      k = 1
      Do 110 j = 1, 20
         Do 100 i = 1, 20, 1
            Qij(i,j) = 0.0001 * Float( IQ( k ) ) 
            k = k + 1
  100       continue
  110    continue
      Return
      End
C
C
      Subroutine dp20( Bits, Qij, n20, MatDef )
      Implicit None
      Integer       i, j, n20
      Real*4        Match, MisMch, Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  Match  / 0.8234  /,   MisMch   /  0.0588 /
      Data  IBits  /  '1.05'  /
      Data  Def(1:40)  / 'Nucleotide PAM 20 transition frequency M' /,
     *      Def(41:60) / 'atrix, Entropy: 1.05' /
C
C
C ***   These Match (9) and MisMatch (-11) values have a relative error
C ***   of 0.0058 with respect to the published values of 1.72 and -2.09
C
C ***   The published values are taken from:  States, D.J., Gish, W., and
C ***   Altschul, S.F.  1991.  Improved Sensitivity of Nucleic Acid
C ***   Database Searches Using Application-Specific Scoring Matrices.
C ***   Methods:  A Companion to Methods in Enzymology, vol. 3,
C ***   pp. 66-70. (August).
C
      MatDef = Def
      Bits = IBits
      Qij(1,1) = Match
      Do 110 i = 2, 4
         Do 100 j = 1, i-1, 1
            Qij(i,j) = MisMch
            Qij(j,i) = MisMch
  100       continue
         Qij(i,i) = Match
  110    continue
      Return
      End
C
C
C
      Subroutine dp50( Bits, Qij, n20, MatDef )
C
      Implicit None
C
      Integer       i, j, n20
      Real*4        Match, MisMch, Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  Match  /  0.6333  /,   MisMch   /  0.1221 /
      Data  IBits  /  '.470'  /
C
      Data  Def(1:40)  / 'Nucleotide PAM 50 transition frequency M' /,
     *      Def(41:60) / 'atrix, Entropy: 0.47' /
C
C ***   These Match (9) and MisMatch (-7) scores have a relative error
C ***   of 0.0021 with respect to the published values of 1.34 and -1.04.
C ***   The scores are log-odds scores for DNA that have diverged to a
C ***   level of PAM 50  (50 point accepted mutations per 100 nucleotides)
C
      MatDef = Def
      Bits = IBits
      Qij(1,1) = Match
      Do 110 i = 2, 4
         Do 100 j = 1, i-1, 1
            Qij(i,j) = MisMch
            Qij(j,i) = MisMch
  100       continue
         Qij(i,i) = Match
  110    continue
      Return
      End
C
      Subroutine dp80( Bits, Qij, n20, MatDef )
C
      Implicit None
C
      Integer       i, j, n20
      Real*4        Match, MisMch, Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  Match  / 0.5061 /,   MisMch   / 0.1644 /
      Data  IBits  /  '0.19'  /
C
      Data  Def(1:40)  / 'Nucleotide PAM 80 transition frequency M' /,
     *      Def(41:60) / 'atrix, Entropy: 0.19' /
C
      MatDef = Def
      Bits = IBits
      Qij(1,1) = Match
      Do 110 i = 2, 4
         Do 100 j = 1, i-1, 1
            Qij(i,j) = MisMch
            Qij(j,i) = MisMch
  100       continue
         Qij(i,i) = Match
  110    continue
      Return
      End
C
      Subroutine dp110( Bits, Qij, n20, MatDef )
C
      Implicit None
C
      Integer       i, j, n20
      Real*4        Match, MisMch, Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  Match  /  0.4211  /,   MisMch   /  0.1926  /
      Data  IBits  /  '0.10'  /
C
      Data  Def(1:40)  / 'Nucleotide PAM 110 transition frequency ' /,
     *      Def(41:60) / 'Matrix, Entropy: .10' /
C
      MatDef = Def
      Bits = IBits
      Qij(1,1) = Match
      Do 110 i = 2, 4, 1
         Do 100 j = 1, i-1, 1
            Qij(i,j) = MisMch
            Qij(j,i) = MisMch
  100       continue
         Qij(i,i) = Match
  110    continue
      Return
      End
C
      Subroutine dp20b( Bits, Qij, n20, MatDef )
C
      Implicit None
C
      Integer       i, j, k, n20
      Integer       Match, Mis( 6 ), XX
      Real*4        Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
C         A
C         C   1                    Subscripting for Array MIS( 6 )
C         G   2     3
C        T/U  4     5     6
C             A     C     G    T/U
C
      Data  Mis  / 0.0371, 0.1008, 0.0371, 0.0371, 0.1008, 0.0371 /
      Data  Match  / 0.8250 /,   XX   /  0.0371   /
      Data  IBits  /  '1.09'  /
C
      Data  Def(1:40)  / '3 to 1 DNA PAM 110 transition frequency ' /,
     *      Def(41:60) / 'Matrix, Entropy: 1.1' /
C      Data  Def(1:40)  / 'DNA PAM 20, 1 transversion / 3 transitio' /,
C     *      Def(41:60) / 'ns, Entropy = 1.09  ' /
C
C
      MatDef = Def
      Bits = IBits
      Qij(1,1) = Match
      k = 0
      Do 110 i = 2, 4
         Do 100 j = 1, i-1, 1
            k = k + 1
            Qij(i,j) = Mis( k )
            Qij(j,i) = Mis( k )
  100       continue
         Qij(i,i) = Match
  110    continue
      Return
      End
C
      Subroutine dp50b( Bits, Qij, n20, MatDef )
C
      Implicit None

C
      Integer       i, j, k, n20
      Real*4        Qij(n20,n20), Match, Mis( 6 ), XX
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  Mis  / 0.0827, 0.1941, 0.0827, 0.0827, 0.1941, 0.0827 /
      Data  Match  / 0.6405  /,   XX   /  0.0827  /
      Data  IBits  /  '0.53'  /
C
      Data  Def(1:40)  / ' 3 to 1 DNA PAM 50 transition frequency ' /,
     *      Def(41:60) / 'Matrix, Entropy: .53' /
C
C
      MatDef = Def
      Bits = IBits
      Qij(1,1) = Match
      k = 0
      Do 110 i = 2, 4, 1
         Do 100 j = 1, i-1, 1
            k = k + 1
            Qij(i,j) = Mis( k )
            Qij(j,i) = Mis( k )
  100       continue
         Qij(i,i) = Match
  110    continue
      Return
      End
C
      Subroutine dp80b( Bits, Qij, n20, MatDef )
C
      Implicit None
C
      Integer       i, j, k, n20
      Real*4        Match, Mis( 6 ), Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  Mis  / 0.1185, 0.2439, 0.1185, 0.1185, 0.2439, 0.1185 /
      Data  Match  /  0.5191  /
      Data  IBits  /  '.283'  /
C
      Data  Def(1:40)  / ' 3 to 1 DNA PAM 80 transition frequency ' /,
     *      Def(41:60) / 'Matrix, Entropy: .28' /
C
C
      MatDef = Def
      Bits = IBits
      Qij(1,1) = Match
      k = 0
      Do 110 i = 2, 4
         Do 100 j = 1, i-1, 1
            k = k + 1
            Qij(i,j) = Mis( k )
            Qij(j,i) = Mis( k )
  100       continue
         Qij(i,i) = Match
  110    continue
      Return
      End
C
C
      Subroutine dp110b( Bits, Qij, n20, MatDef )
C
      Implicit None
C
      Integer       i, j, k, n20
      Real*4        Match, Mis( 6 ), Qij(n20,n20)
      Character*4   Bits, IBits
      Character*60  MatDef, Def
C
      Data  Mis  / 0.1467, 0.2685, 0.1467, 0.1467, 0.2685, 0.1467 /
      Data  Match  /  0.4381  /
      Data  IBits  /  '.157'  /
C
      Data  Def(1:40)  / '3 to 1 DNA PAM 110 transition frequency ' /,
     *      Def(41:60) / 'Matrix, Entropy: .16' /
C
C
      MatDef = Def
      Bits = IBits
      Qij(1,1) = Match
      k = 0
      Do 110 i = 2, 4
         Do 100 j = 1, i-1, 1
            k = k + 1
            Qij(i,j) = Mis( k )
            Qij(j,i) = Mis( k )
  100       continue
         Qij(i,i) = Match
  110    continue
      Return
      End
C
      SUBROUTINE PAKNAM ( NAME, LENGTH, BLANKS )
C
C ****************************************************************************
C **   Copyright  --  Pittsburgh Supercomputing Center  --  December 1988   **
C ********************  Written by Hugh B. Nicholas  *************************
C
C ***   Subroutine PAKNAM (PAcK NAMe) removes leading spaces from a string
C ***   of text stored in a character variable.
C
      IMPLICIT NONE
C
C ***   Passed variables
C
      INTEGER        LENGTH      !   position of last non-space character
      CHARACTER*(*)  NAME        !   text string
      LOGICAL        BLANKS      !   set .TRUE. if text string is all spaces
C
C ***   Local Variables
C
      INTEGER        I, J, L
C
      BLANKS = .FALSE.
      DO 100 L = 1, LENGTH
         IF( NAME(L:L) .EQ. ' ' )    GO TO 100
            I = L
            GO TO 150
  100    CONTINUE
      BLANKS = .TRUE.
      RETURN                            !   string was all spaces
C
  150 IF( I .GT. 1 )    THEN
         J = 0
         DO 200 L = I, LENGTH
            J = J + 1
            NAME(J:J) = NAME(L:L)
  200       NAME(L:L) = ' '
      ELSE
      END IF
C
      RETURN
      END
C
C
      INTEGER FUNCTION POSN( SIZE, STRING, TABLE )
C
C ****************************************************************************
C **   Copyright  --  Pittsburgh Supercomputing Center  --  December 1988   **
C ********************   Coded by Hugh B. Nicholas   *************************
C
C
C ***   The functions, STDPSN is a binary search of the table STDPOS .
C ***   search algorithm is an implementation of Knuth's algorithm U
C ***   (vol 3., ch. 6.2.1, p 411; see vol 1., ch. 1.2.4, p 37 for
C ***   details of the notation).
C
C ***   INDEX = The current position in the search table.
C ***   RATE  = The change in INDEX if a match is not found on the current
C ***         cycle.  If RATE becomes zero the value passed to the
C ***         function is not in the search table.
C
C
      INTEGER         RATE, INDEX, SIZE
      CHARACTER*(*)   STRING, TABLE( 0 : SIZE )
C
C
      RATE = SIZE  / 2
      INDEX = ( SIZE + 1 ) / 2
  200 IF( LLT( STRING, TABLE( INDEX ) ) )    THEN
         INDEX = INDEX - ( ( RATE + 1 ) / 2 )
      ELSE IF( LGT( STRING, TABLE( INDEX ) ) )    THEN
         INDEX = INDEX + ( ( RATE + 1 ) / 2 )
      ELSE IF( STRING .EQ. TABLE( INDEX ) )    THEN
         POSN = INDEX
         RETURN
      END IF
      IF( RATE .EQ. 0 )    THEN
         POSN = 0
         RETURN
      ELSE
      END IF
      RATE = RATE / 2
      GO TO 200
C
      END
C
C
      INTEGER FUNCTION LASTCH ( STRING, MAX )
C
C ****************************************************************************
C **   Copyright  --  Pittsburgh Supercomputing Center  --  December 1988   **
C ********************  Written by Hugh B. Nicholas  *************************
C
C
C ***   FUNCTION LASTCH finds the last non-space character in a string
C ***   stored in a character variable.
C
      IMPLICIT NONE
C
      INTEGER         MAX, I
      CHARACTER*(*)   STRING
C
      DO 100 I = MAX, 1, -1
         IF( STRING(I:I) .EQ. ' ' )    GO TO 100
            LASTCH = I
            RETURN
  100       CONTINUE
C
      LASTCH = 0
      RETURN               !   the string was all blanks
C
      END
C
C
      SUBROUTINE ALLCAP( NAME, LENGTH )
C
C  ****************************************************************************
C  **      Copyright   Pittsburgh Supercomputing Center  --  May 1988        **
C  ********************   written by  Hugh B. Nicholas   **********************
C
C ***   Subroutine ALLCAP converts the alphameric characters in the first
C ***   LENGTH characters of the character variable NAME to upper case.
C
      IMPLICIT NONE
C
      INTEGER        LENGTH, L
      CHARACTER*(*)  NAME
C
      DO 100 L = 1, LENGTH
         IF( LGE( NAME(L:L), 'a' )   .AND.   LLE( NAME(L:L), 'z' ) )
     *       NAME(L:L) = CHAR( ICHAR( NAME(L:L) ) - 32 )
  100    CONTINUE
C
      RETURN
      END
C
C
      INTEGER FUNCTION STYLE( NUMBER, LENGTH )
C
C ****************************************************************************
C **   Copyright  --  Pittsburgh Supercomputing Center  --  December 1988   **
C ********************  Written by Hugh B. Nicholas  *************************
C
C
C ***   Function style determines the type of FORMAT descriptor is necessary
C ***   to read the numeric value stored as ASCII characters in the variable
C ***   NUMBER.  Error checking is fairly rigorous for both the content and
C ***   format of elements to be read as numeric data.  An ANSI standard
C ***   FORTRAN fixed field numeric format read allows imbedded blanks
C ***   around signs (+,-) and around exponent indicators (E,D).  For
C ***   example,  - 12.9735 E - 01   will be read correctly with an E16.x
C ***   format as  -0.129735E+01.  Because this program used spaces (blanks)
C ***   as delimiters to allow free format, the imbedded blanks will cause
C ***   the number to be divided into five tokens.  Also entries such as
C ***   .E  or  E   are read as zero in a fixed format field while this
C ***   program will not translate a field as numeric unless it contains at
C ***   least one digit, i.e.,  0.E   .0E or   0E    will all be read as
C ***   zero by this program.
C
      INTEGER         LENGTH
      CHARACTER*31    NUMBER
C
      INTEGER   NLEGLP, TYPES
C
      PARAMETER  (  NLEGLP = 17,  TYPES = 31 )
C
      INTEGER         I, DIGIT, D, E, F, BAD, DCOUNT, ECOUNT, SCOUNT,
     *                FCOUNT, FPOS, EPOS, SPOS1, SPOS2, DIGIT1
      INTEGER         FORM( 0 : TYPES )
      CHARACTER*1     LEGALC( 0 : NLEGLP )
      INTEGER         POSN                                  !   function name
C
C       STYLE           Indicates numeric format of type:
C      returns
C
C        -1                  An illegal NUMERIC character is present in
C                               the field or the format is incorrect, read
C                               as alphameric.
C         1                  Integer
C         2                  Real - has decimal point but no exponent
C         3                  Exponential - has an E exponent
C         4                  Double Precision Exponential - has a D exponent
C
C ***   LEGALC contains the ASCII characters which may be legally used
C ***   to represent a number for FORTRAN-77 input and output.  The first
C ***   element, ' ', is a sentinal needed by the binary search (POSN)
C ***   subroutine.  The punctuation characters are in ASCII collating
C ***   sequence order.  The letter 'Q' used by VAX fortran as an
C ***   exponent indicator for REAL*16 numbers is not included.
C
      DATA LEGALC / ' ', '+', '-', '.', '0', '1', '2', '3', '4', '5',
     +              '6', '7', '8', '9', 'D', 'E', 'd', 'e'            /
C
      DATA FORM  /  -1, 1, 2, 2, -1, 3, -1, 3, -1, 4, -1, 4, 20 * -1 /
C
      BAD = 0
      DIGIT = 0
      D = 0
      E = 0
      F = 0
      FCOUNT = 0
      DCOUNT = 0
      ECOUNT = 0
      SCOUNT = 0
      FPOS = 0
      EPOS = 0
      SPOS1 = 0
      SPOS2 = 0
      DIGIT1 = 0
      DO 100 I = 1, LENGTH
         CHTYPE = POSN( NLEGLP, NUMBER(I:I), LEGALC )
         IF( CHTYPE .GE. 4   .AND.   CHTYPE .LE. 13 )    THEN
            IF( DIGIT .EQ. 0 )    DIGIT1 = I
            DIGIT = 1
         ELSE IF( CHTYPE .EQ. 3 )    THEN
            F = 2
            FCOUNT = FCOUNT + 1
            FPOS = I
         ELSE IF( CHTYPE .EQ. 15   .OR.   CHTYPE .EQ. 17 )    THEN
            E = 4
            ECOUNT = ECOUNT + 1
            EPOS = I
         ELSE IF( CHTYPE .EQ. 14   .OR.   CHTYPE .EQ. 16 )    THEN
            D = 8
            DCOUNT = DCOUNT + 1
            EPOS = I
         ELSE IF( CHTYPE .EQ. 1   .OR.   CHTYPE .EQ. 2 )    THEN
            SCOUNT = SCOUNT + 1
            IF( SCOUNT .EQ. 1 )    THEN
               SPOS1 = I
            ELSE IF( SCOUNT .EQ. 2 )    THEN
               SPOS2 = I
            ELSE
            END IF
         ELSE
            BAD = 16
         END IF
  100    CONTINUE
C
C ***   The style indicators have been set now get the overall style or
C ***   format type.
C
      STYLE = FORM( DIGIT + D + E + F + BAD )
C
C ***   If the style check indicates that the data is in a numeric format,
C ***   do a more complete format check.
C
      IF( STYLE .EQ. 1 )    THEN
         IF( SCOUNT .GT. 1   .OR.   SPOS1 .GT. DIGIT1 )    STYLE = -1
      ELSE IF( STYLE .EQ. 2 )    THEN
         IF( SCOUNT .GT. 1        .OR.
     *       FCOUNT .GT. 1        .OR.
     *       SPOS1  .GT. DIGIT1   .OR.
     *       SPOS1  .GT. FPOS )           STYLE = -1
      ELSE IF( STYLE .EQ. 3   .OR.   STYLE .EQ. 4 )    THEN
         IF( SCOUNT .GT. 2     .OR.
     *       FCOUNT .GT. 1     .OR.
     *       FPOS   .GT. EPOS )           STYLE = -1
         IF( SCOUNT .EQ. 0 )    THEN
         ELSE IF( SCOUNT .EQ. 1 )    THEN
            IF( SPOS1 .GT. DIGIT1 .AND. SPOS1 .NE. EPOS+1 )  STYLE = -1
         ELSE IF( SCOUNT .EQ. 2 )    THEN
            IF( SPOS1 .GT. DIGIT1  .OR.
     *          SPOS1 .GT. EPOS    .OR.
     *          SPOS2 .NE. EPOS+1 )         STYLE = -1
         ELSE
            STYLE = -1
         END IF
      ELSE
      END IF
C
      RETURN
      END
C
      Subroutine AveVar( Data, N, Ave, Var )
C
      Implicit None
C
      Integer     j, N
      Real*4      Data( N ), Ave, Var
C
      Ave = 0.0
      Var = 0.0
C
      Do 100 j = 1, N
         Ave = Ave + Data( j )
  100    continue
      Ave = Ave / Float( N )
C
      Do 200 j = 1, N
         Var = Var + ( Data( j ) - Ave ) ** 2
  200    continue
      Var = Var / Float( N - 1 )
      Var = Sqrt( Var )
C
      Return
      End
C
C
      Subroutine hpsort( N, RA )
      Integer    N
      Real*4     RA(n)
      Integer    i, j, l, ir
      Real*4     rra
C
      If( N .lt. 2 )    Return
      l = N / 2 + 1
      ir = N
  100 If( l .gt. 1 )    Then
         l = l - 1
         rra = ra(l)
      Else
         rra = ra(ir)
         ra(ir) = ra(1)
         ir = ir - 1
         If( ir .eq. 1 )    Then
            ra(1) = rra
            Return
         Else
         EndIf
      EndIf
      i = l
      j = l + l
  200 If( j .le. ir )    Then
         If( j .lt. ir )    Then
            If( ra(j) .lt. ra(j+1) )    j = j+1
         Else
         EndIf
         If( rra .lt. ra(j) )    Then
            ra(i) = ra(j)
            i = j
            j = j + j
         Else
            j = ir + 1
         EndIf
         GoTo 200
      EndIf
      ra(i) = rra
      GoTo 100
C
      End
C
C
      Subroutine RmGaps( MaxSqP, TotSqP, NSP, NS, MxSymP, AlnLen, Gap,
     *                   LNoGap, LsCum, Seq, GrupId, NoGaps, GapMap,
     1                   NSqUsed, FrGap, BackFr, ASize, UnKnown, GComp,
     3                   OutFl, MaxVol, MinVol, AveVol, SigVol, Consen,
     4                   SqPrnt, Termnl, KeyBrd, GrupFl, SqName )
C
      Implicit None
C
C ***   Passed variable declarations
C
      Integer       MaxSqP, TotSqP, NSP, NS, MxSymP, AlnLen, Gap,
     *              LNoGap, NSqUsed, OutFl, ASize, UnKnown, Termnl,
     *              KeyBrd, GrupFl
      Integer       LsCum( NSP ), Seq( TotSqP ), NoGaps( MaxSqP ),
     *              GapMap( MaxSqP ), GComp( 0:MxSymP ), GrupId( NSP )
      Real*4        FrGap( MaxSqP ), BackFr( MxSymP ), MaxVol( MaxSqP ),
     *              MinVol( MaxSqP ), AveVol( MaxSqP ), SigVol(MaxSqP)
      Character*1   SqPrnt(25), Consen( MaxSqP )
      Character*15  SqName( NSP )
C
C ***   Local Variable declarations
C
C ***          LocNsP must be at least as large as Global Parameter NSP
C ***          LMxSqP must be at least as large as Global Parameter MaxSqP
C
      Integer     LocNsP, LMxSqP
      Parameter ( LocNsP = 600,  LMxSqP = 1000 )
C
      Integer     l, n, aa, tc(0:25), nb, NBlock, i, ib, ie, aaknt,
     *            maxaa, NKeep, KeepIt( LMXSqP ), MaxGap
      Real*4      FUsed
      Real*4      StdVol( 23 ), tvol( LocNsP )
      Save        StdVol
C
C       Amino Acid/ 'a',  'r',   'n',   'd',   'c',   'q',   'e',
      Data StdVol / 90.0, 195.0, 126.0, 118.0, 113.0, 150.0, 142.0,
C
C                   'g',  'h',   'i',   'l',   'k',   'm',   'f',
     *              64.0, 159.0, 164.0, 164.0, 170.0, 167.0, 193.0,
C
C                   'p',   's',  't',   'w',   'y',   'v',   'b',
     1              124.0, 95.0, 121.0, 231.0, 197.0, 139.0, 122.0,
C
C                   'z',    'x',
     2              146.0,  0.0   /
C
    1 Format(/' Mapping of positions used in the analysis (alignment',
     *        ' columns without Gaps)',
     1       /' to the original complete alignment.' )
    2 Format(/' WithOutGaps: ', 12(I5:) )
    3 Format( ' OriginalPos: ', 12(I5:) )
    4 Format( /' Please enter the maximum number of gap characters',
     *         ' allowed in the',
     1        /' alignment positions of the trimmed alignment to be',
     2         ' written.',
     3        /' A negative number will suppress writing the trimmed',
     4         ' alignment.', / )
    5 Format('>', A15 )
    6 Format('   ', 75A1: )
C
C
      NKeep = 0
      Write( Termnl, 4 )
      Read( KeyBrd, * )  MaxGap
      If( MaxGap .ge. 0 )    Then
         Open( Unit = GrupFl, file = 'trimmed.aa', status = 'new',
     *         Form = 'formatted', Access = 'Sequential')
c     *         CarriageControl = 'List' )
      EndIf
      LNoGap = 0
      FUsed = Float( NSqUsed )
      Do 120 n = 1, MxSymP, 1
         GComp(n) = 0
  120    continue
      StdVol(UnKnown) = 0.0
      Do 140 n = 1, ASize, 1
         StdVol(UnKnown) = StdVol(UnKnown) + BackFr(n) * StdVol(n)
  140    continue
C
      Do 400 l = 1, AlnLen, 1
         aaknt = 0
         MaxVol(l) = -10.0
         MinVol(l) = 500.0
         Do 200 n = 1, 25, 1
            tc(n) = 0
  200       continue
         Do 300 n = 1, NS, 1
            If( GrupId(n) .ge. 1 )    Then
               aa = Seq( LsCum(n) + l )
               tc( aa ) = tc( aa ) + 1
               If( aa .ne. gap )    Then
                  aaknt = aaknt + 1
                  tvol(aaknt) = StdVol(aa)
                  If( tvol(aaknt) .gt. MaxVol(l)) MaxVol(l)=tvol(aaknt)
                  If( tvol(aaknt) .lt. MinVol(l)) MinVol(l)=tvol(aaknt)
               Else
               EndIf
            Else
            EndIf
  300       continue
         If( tc(gap) .eq. 0 )    Then
            Do 310 n = 1, 25, 1
               GComp(n) = GComp(n) + tc(n)
  310          continue
            LNoGap = LNoGap + 1
            NoGaps(LNoGap) = l
            GapMap(l) = LNoGap
            FrGap(l) = 0.0
         Else
            GapMap(l) = -tc(gap)
            FrGap(l) = Float( tc(gap) ) / FUsed
            n = 0
            maxaa = 0
            Do 330 aa = 1, ASize, 1
               If( tc(aa) .gt. maxaa )    Then
                  n = aa
                  maxaa = tc(aa)
               Else
               EndIf
  330          continue
            Consen(l) = SqPrnt(n)
         EndIf
C
         If( aaknt .gt. 2 )    Then
            Call AveVar( tvol, aaknt, AveVol(l), SigVol(l) )
         Else
            AveVol(l) = ( MaxVol(l) + MinVol(l) ) / 2.0
            SigVol(l) = -1.0
         EndIf
C
         If( tc(gap) .le. MaxGap   .and.   MaxGap .ge. 0 )    Then
            NKeep = NKeep + 1
            KeepIt( NKeep ) = l
         Else
         EndIf
  400    continue
C
      NBlock = ( LNoGap + 11 ) / 12
      Write( OutFl, 1 )
      ib = -11
      Do 500 nb = 1, NBlock, 1
         ib = ib + 12
         ie = ib + 11
         If( ie .gt. LNoGap )    ie = LNoGap
         Write( OutFl, 2 )  ( i, i = ib, ie, 1 )
         Write( OutFl, 3 )  ( NoGaps(i), i = ib, ie, 1 )
  500    continue
C
      If( MaxGap .ge. 0 )    Then
         Do 700 n = 1, NS, 1
            Write( GrupFl, 5 )  SqName(n)
CAJR2020
C            write(6,*) 'At 6292',lsCum(n),NKeep
C            write(6,*) 'at 6293',(Keepit(l),l=1,Nkeep,1)
C            write(6,*) 'at 6294',(Seq(LsCum(n)+Keepit(l)),l=1,Nkeep,1)
            Write( GrupFl, 6 )  ( SqPrnt( Seq( LsCum(n) + KeepIt(l) )),
     *                           l = 1, NKeep, 1 )
  700    continue
         Close( Unit = GrupFl, Status = 'Keep' )
      Else
      EndIf
C
      Return
      End
C
C
      Subroutine TwoOut( Termnl, GrupFl, NSPSqP, NSP, WKnt, BKnt,
     *                   PFLen, PFile, NSqUsed, SeqUsed, SqName, WAve,
     *                   BAve, Wstd, Bstd, Within, Between, Pair )
C
      Integer         Termnl, GrupFl, NSPSqP, NSP, WKnt, BKnt, PFLen,
     *                NSqUsed
      Integer         SeqUsed( NSP )
      Real*4          WAve, BAve, Wstd, Bstd
      Real*4          Within(NSPSQP), Between( NSPSQP ),
     *                Pair( NSP, NSP )
      Character*15    SqName( NSP )
      Character*24    PFile
C
C ***   Local varialbles -  Warning:  LNPSSQ should be the same as
C       the global parameter NSPSQP = (( NSP * ( NSP + 1 )) / 2 )
C
      Integer        LNSP, LNSPSQ
      Parameter   (  LNSP = 600,  LNSPSQ = ((LNSP*(LNSP+1))/2 )  )
C
      Integer         n, na, NBlock, i, ib, ie
      Real*4          FWKnt, FBKnt
      Real*4          WFrac( LNSPSQ ), BFrac( LNSPSQ )
      Logical         Skip
      Character*34    CVFile
C
   35 Format(//' The local parameter "LNSP" in subroutine TwoOut is',
     *         ' too small.  It must be',
     1        /' set so that ( LNSP * ( LNSP + 1 )) / 2  is at least:',
     2           I12,
     3        /' Because of this the Within and Between distance',
     4         ' distribution files will',
     5        /' not be produced by this run.', // )
   36 Format(/' The Within Groups Distribution:',
     *       /'                        Members:', I10,
     1       /'                        Average:', F10.4,
     2       /'             Standard Deviation:', F10.4 )
   37 Format(/' The Between Groups Distribution:',
     *       /'                         Members:', I10,
     1       /'                         Average:', F10.4,
     2       /'              Standard Deviation:', F10.4 )
   38 Format(' Count,  Within,  Fraction')
   39 Format(' ', I4, ',', F8.3, ',', F7.4 )
   40 Format(' Count,  Between,  Fraction')
   41 Format(' ', I10 )
   42 Format( A10, ' ', 10( F10.3: ) )
   43 Format( 11X, 10( F10.3: ) )
C
      Write( Termnl, 36 )  WKnt, WAve, Wstd
      Write( Termnl, 37 )  BKnt, BAve, Bstd
      Skip = .False.
      If( WKnt .gt. LNSPSQ )    Then
         Write( Termnl, 35 )  WKnt
         Skip = .True.
      Else
      EndIf
      If( BKnt .gt. LNSPSQ )    Then
         Write( Termnl, 35 )  BKnt
         Skip = .True.
      Else
      EndIf
      If( Skip )    GoTo 300
C
C
      Open( Unit = GrupFl, File = 'WithIn.Plt', Status = 'new',
     *      Access = 'Sequential', form = 'Formatted')
c     *      CarriageControl = 'List' )
      FWKnt = 1.0 / Float( WKnt )
      Call HpSort( WKnt, WithIn )
      Write( GrupFl, 38 )
      Do 260 n = 1, WKnt, 1
         WFrac( n ) = Float( n ) * FWKnt
         Write( GrupFl, 39 )  n, Within(n), WFrac(n)
  260    continue
      Close( Unit = GrupFl, Status = 'keep' )
C
      Open( Unit = GrupFl, File = 'Between.Plt', Status = 'new',
     *      Access = 'Sequential', form = 'Formatted')
c     *      CarriageControl = 'List' )
      FBKnt = 1.0 / Float( BKnt )
      Call HpSort( BKnt, Between )
      Write( GrupFl, 40 )
      Do 280 n = 1, BKnt, 1
         BFrac( n ) = Float( n ) * FBKnt
         Write( GrupFl, 39 )  n, Between(n), BFrac(n)
  280    continue
      Close( Unit = GrupFl, Status = 'keep' )
C
C
C ***   Write the Symetrical Pairwise Entropy Distances into a file
C
  300 CVFile = ' '
      CVFile(1:PfLen) = PFile(1:PfLen)
      CVFile(PfLen+1:PfLen+4) = '.dis'
      Open( Unit = GrupFl, File = CVFile, Status = 'new',
     *      Access = 'Sequential', form = 'Formatted')
c     *      CarriageControl = 'List' )
      Write( GrupFl, 41 )  NSqUsed
      NBlock = ( NSqUsed + 9 ) / 10
      Do 500 na = 1, NSqUsed, 1
         n = SeqUsed( na )
         ib = 1
         ie = 10
         If( ie .gt. NSqUsed )    ie = NSqUsed
C         Write( GrupFl, * )  ' na = ', na
C         Write( GrupFl, * )  ' n = ', n
C         Write( GrupFl, * )  ' SqName(n) = ', SqName(n)
C         Write( GrupFl, * )  ' ib = ', ib
C         Write( GrupFl, * )  ' SeqUsed(ib) = ', SeqUsed(ib)
C         Write( GrupFl, * )  ' Pair(n,SeqUsed(ib)) = ',
C     *                         Pair(n,SeqUsed(ib))
C         Write( GrupFl, * )  ' ie = ', ie
C         Write( GrupFl, * )  ' SeqUsed(ie) = ', SeqUsed(ie)
C         Write( GrupFl, * )  ' Pair(n,SeqUsed(ie)) = ',
C     *                         Pair(n,SeqUsed(ie))
         Write( GrupFl, 42 )  SqName(n)(1:10), ( Pair( n, SeqUsed(i) ),
     *                          i = ib, ie, 1 )
         Do 450 nb = 2, NBlock, 1
            ib = ib + 10
            ie = ib + 9
            If( ie .gt. NSqUsed )    ie = NSqUsed
            Write( GrupFl, 43 )  ( Pair(n,SeqUsed(i)), i = ib, ie, 1 )
  450       continue
  500    continue
C
      Return
      End
C
      Subroutine IAsnGp( Termnl, Title, GrpNam, GrpKnt, NGroup, InvGrp,
     *                   SqName, Fscore, GrpAln, GrupId, NS, NSP,
     1                   MxGrpP, MaxSqP, NSqUsed, SeqUsed, KeyBrd,
     2                   Profil, ProfileT, GZEnt, MaxGZ, ZLimit,
     3                   TotSqP, LNoGap, NoGaps, LsCum, Seq )
C
C ***   Passed Variable declarations
C
      Integer       NSP, MxGrpP, NS, Termnl, KeyBrd, GrpKnt, NSqUsed,
     *              TotSqP, LNoGap
      Integer       Ngroup( MxGrpP ), InvGrp( NSP, MxGrpP ),
     *              GrupId(NSP), SeqUsed( NSP ), LsCum( NSP ),
     *              Seq( TotSqP ), NoGaps( MaxSqP )
      Real*4        FScore( NSP ), GrpAln( NSP, 0:MxGrpP ),
     *              Profil( 25, MaxSqP, 2, MxGrpP ), MaxGZ( MaxSqP ),
     *              GZEnt( MaxSqP, MxGrpP ), ProfileT( 23, MaxSqP ),
     *              ZLimit
      Character*15  SqName( NSP )
      Character*30  GrpNam( MxGrpP )
      Character*70  Title
C
C ***   Local variable declarations
C ***   Note:  Array NewKnt( MGpLoc ) should be dimensioned to be at
C ***   least as large as the Parameter  MxGrpP  in the main program.
C
      Integer      MGpLoc
      Parameter  ( MGpLoc = 25 )
C
      Integer      DoGrup, Nblock, l, la, ib, ie, n, sn, HiGrup,
     *             NewGrp, NewMem, aatype
      Integer      NewKnt(0:MGpLoc)
      Real*4       HiScor, FamScr, GrpScr( MGpLoc )
      Logical      UseAll
C
    1 Format(/' Sequence: ', A15, 10X, 'Alignment Score with Family =',
     *           F11.3,
     1       /' High Scoring Group is: ', A30, ',   Score =', F11.3 )
    2 Format(/' Group Number, Name, and Score', 11X,
     *        ' Group Number, Name, and Score' )
    3 Format( I4,'. ', A20, ' =', F10.2:, I7, '. ', A20, ' =', F10.2 ) 
    4 Format(/' Enter the number of the group to assign sequence ',A15,
     *        ' to a',
     1       /' specific group.  Enter zero (0) to go to the next',
     2        ' sequence without',
     3       /' assigning the current sequence to a group, or enter',
     4        ' minus 1 (-1) to',
     5       /' quit assigning sequences to groups and return to the',
     6        ' main program.' )
    5 Format(/' Sequence ', A15, ' is now in group: ', A30 )
    6 Format(/' Sequence ', A15, ' is not assigned to any group.' )
    7 Format(/' ', I5, ' is not a recognized choice, please re-enter.')
    8 Format(/' All of the sequences initially in the To Be Assigned',
     *        ' Group have been',
     1       /' assigned to one of original defined groups' )
    9 Format(/ I5, ' sequences are still in the To Be Assigned Group.')
C
      If( ZLimit .lt. 0 )    Then
         UseAll = .True.
      Else
         UseAll = .False.
      EndIf
      NewMem = 0
      Do 100 n = 0, MGpLoc, 1
         NewKnt( n ) = 0
  100    continue
C
C ***   Write the report of data for assigning previously unclassified
C ***   sequences to specific groups to the terminal screen one
C ***   sequence at a time.
C ***   1000 loop over all sequences in the highest numbered formal group.
C ***   This group contains the sequences designated by the user to be
C ***   assigned to one of the user defined groups of sequences.
C
      Do 1000 n = 1, NGroup( GrpKnt ), 1
         sn = InvGrp(n,GrpKnt)
         If( UseAll )    Then
            FamScr = FScore( sn )
            Do 200 DoGrup = 1, GrpKnt - 1, 1
               GrpScr( DoGrup ) = GrpAln( sn, DoGrup )
  200          continue
         ElseIf( .not. UseAll )    Then
            FamScr = 0.0
            Do 250 DoGrup = 1, GrpKnt - 1, 1
               GrpScr( DoGrup ) = 0.0
  250          continue
            Do 400 la = 1, LNoGap, 1
               l = NoGaps(la)
               If( MaxGZ(l) .ge. ZLimit )    Then
                  aatype = Seq( LsCum(sn)+l )
                  FamScr = FamScr + ProfileT( aatype, l )
                  Do 300 DoGrup = 1, GrpKnt - 1, 1
                     GrpScr( DoGrup ) = GrpScr(DoGrup)
     *                                + Profil( aatype, l, 1, DoGrup )
  300                continue
               Else
               EndIf
  400          continue
         EndIf
C
C ***   Find the group whose profile gives the highest score with the
C ***   current sequence.
C
         HiGrup = 1
         HiScor = GrpScr( 1 )
         Do 700 DoGrup = 2, GrpKnt - 1, 1
            If( GrpScr( DoGrup )  .gt.  HiScor )    Then
               HiScor = GrpScr( DoGrup )
               HiGrup = DoGrup
            Else
            EndIf
  700       continue
         Write( Termnl, 1 )  SqName(sn), FamScr, GrpNam( HiGrup ),
     *                       HiScor
C
C ***   List the scores with profiles for all groups on the screen
C
         Write( Termnl, 2 )
         NBlock = GrpKnt / 2
         ib = -1
         Do 750 nb = 1, NBlock, 1
            ib = ib + 2
            ie = ib + 1
            If( ie .gt. GrpKnt-1 )    ie = GrpKnt-1
CAJR - commented out to get around compiler error
C            Write( Termnl, 3 ) (( DoGrup, GrpNam(DoGrup)(1:20),
C     *                          GrpScr(DoGrup) ), DoGrup = ib, ie, 1)
  750       continue
C
C ***   Accept the assignment of the current sequence to a group by
C ***   the user or skip to the next sequence to be assigned.
C
         Write( Termnl, 4 )  SqName( sn )
  800    Read( KeyBrd, * )   NewGrp
         If( NewGrp .ge. 1   .and.   NewGrp .le. (GrpKnt-1) )   Then
            Write( Termnl, 5 )  SqName( sn ), GrpNam( NewGrp )
            GrupID( sn ) = NewGrp
            NewKnt( NewGrp ) = NewKnt(NewGrp) + 1
            NewKnt( GrpKnt ) = NewKnt(GrpKnt) - 1
            NewMem = NewMem + 1
         ElseIf( NewGrp .eq. 0 )    Then
            Write( Termnl, 6 )  SqName( sn )
         ElseIf( NewGrp .le. -1 )   Then
            Write( Termnl, 6 )  SqName( sn )
            GoTo 1100
         Else
            Write( Termnl, 7 )  NewGrp
            GoTo 800
         EndIf
C
 1000    continue
C
 1100 If( NewMem .gt. 0 )    Then
         Do 1200 DoGrup = 1, GrpKnt, 1
            Do 1150 n = 1, NGroup(DoGrup), 1
               InvGrp(n,DoGrup) = 0
 1150          continue
            NGroup( DoGrup ) = NGroup(DoGrup) + NewKnt(DoGrup)
 1200       continue
         Do 1300 n = 0, MGpLoc, 1
            NewKnt( n ) = 0
 1300       continue
         Do 1400 n = 1, NS, 1
            DoGrup = GrupID(n)
            NewKnt(DoGrup) = NewKnt(DoGrup) + 1
            If( DoGrup .ge. 1 )    Then
               InvGrp( NewKnt(DoGrup), DoGrup ) = n
               NSqUsed = NSqUsed + 1
               SeqUsed( NSqUsed ) = n
            Else
            EndIf
 1400       continue
         If( NewKnt(GrpKnt) .gt. 0 )    Then
            Write( Termnl, 9 )  NewKnt(GrpKnt)
         Else
            Write( Termnl, 8 )
         EndIF
      Else
      EndIf
C
      Return
      End
C
C
      Subroutine WrGpFl( Termnl, KeyBrd, GrupFl, NSP, MxGrpP, SqName,
     *                   GrpKnt, GrpNam, Title, PFile, InvGrp, NGroup )
C
C ***   Passed Variable declarations
C
      Integer       NSP, MxGrpP, GrupFl, GrpKnt, Termnl, KeyBrd
      Integer       Ngroup( MxGrpP ), InvGrp( NSP, MxGrpP )
      Character*15  SqName( NSP )
      Character*24  PFile
      Character*30  GrpNam( MxGrpP )
      Character*70  Title
C
C ***   Local varaible declarations
C
      Integer       DoGrup, NLine, l, i, ib, ie, n, sn
      Character*80  GpFile, SqLine
      Logical       empty
C
    1 Format(//' Enter the name of a file to hold the new group',
     *         ' assignments.' )
    2 Format( A80 )
    3 Format('Problem  ', A24 )
    4 Format('Title  ', A70 )
    5 Format('Group  ', A30 )
    6 Format(' ', 4( ' ', A15: ) )
    7 Format(//' The new group assignments will be in:  NewGroup.dat')
C
      Write( Termnl, 1 )
      Read( KeyBrd, 2 )  GpFile
      Call PakNam( GpFile, 80, empty )
      If( empty )    Then
         GpFile(1:12) = 'NewGroup.dat'
         Write( Termnl, 7 )
      Else
      EndIf
      Open( Unit = GrupFl, File = GpFile, status = 'new',
     *      Access = 'Sequential', Form = 'Formatted')
c     *      CarriageControl = 'List' )
      Write( GrupFl, 3 )  PFile
      Write( GrupFl, 4 )  Title
      Do 500 n = 1, GrpKnt, 1
         Write( GrupFl, 5 )  GrpNam( n )
         NLine = ( NGroup(n) + 3 ) / 4
         ib = -3
         Do 400 l = 1, NLine, 1
            SqLine = ' '
c            Line = ' '
            ib = ib + 4
            ie = ib + 3
            If( ie .gt. NGroup(n) )    ie = NGroup(n)
            Write( SqLine, 6 )  ( SqName( InvGrp( i, n )), i = ib,ie,1)
            If( l .eq. NLine )    SqLine(71:72) = '//'
            Write( GrupFl, 2 )  SqLine
  400       continue
  500    continue
C
      Return
      End
C
C
      Subroutine ZScore( OutFl, MaxSqP, MxGrpP, NSP, TotSqP,
     *                   Seq, LsCum, AlnLen, LGpNam, Entropy,
     1                   Cross, GrpNam, SqPrnt, LNoGap, GrpKnt,
     2                   InvGrp, SqName, GapMap, GZav, GZstd, FZav,
     3                   FZstd, MaxGZ, FZEnt, GZEnt, NotId )
      Implicit None
C
C ***   Passed Variable Declarations
C
      Integer        OutFl, MaxSqP, MxGrpP, NSP, TotSqP, AlnLen,
     *               LNoGap, GrpKnt, NotId
      Integer        Seq( TotSqP ), LsCum( NSP ), LGpNam( MxGrpP ),
     *               GapMap( MaxSqP ), InvGrp( NSP, MxGrpP )
      Real*4         GZav, GZstd, FZav, FZstd
      Real*4         Entropy( MaxSqP ), Cross(MaxSqP, 2, MxGrpP ),
     *               MaxGZ(MaxSqP), FZEnt(MaxSqP), GZEnt(MaxSqP,MxGrpP)
      Character*1    SqPrnt( 25 )
      Character*15   SqName( NSP )
      Character*30   GrpNam( MxGrpP )
C
C ***   Local Variable Declarations
C
C ***   Note:  The local parameter LMxSqP should be at least as large
C              as the Global parameter MaxSqP.
C
      Integer        LMxSqP
      Parameter  (   LMxSqP = 1000 )
C
      Integer        ng, sn, l, lm
      Integer        GRank( LMxSqP ), FRank( LMxSqP )
      Real*4         ge( LMxSqP ), gp( LMxSqP ),
     *               fe( LMxSqP ), fp( LMxSqP )
C
    5 Format( '//           Group: ', A30,
     *       /'// Master Sequence: ', A15,
     1       /'// Average Family Entropy =', F8.3, ',  Standard',
     2        ' Deviation =', F8.3,
     3       /'//  Average Group Entropy =', F8.3, ',  Standard',
     4        ' Deviation =', F8.3 )
C    7 Format(' Ntop = ', I5, ', Top10, Top20, Top30 =', 3F8.3 )
C    8 Format('   L    Fe    fz    fp      Ge    gz    gp seq  GapMap')
C    9 Format( I5, 6F8.3, 2X, A1, I8 )
C
      Do 200 l = 1, AlnLen, 1
         fe(l) = Entropy(l)
         MaxGZ( l ) = 0.0
  200    continue
      Call Rindex( AlnLen, fe, FRank )
      Do 220 l = 1, LNoGap, 1
         lm = AlnLen + 1 - l
         fp(l) = fe(FRank(lm))
  220    continue
      Call AveVar( fp, LNoGap, FZav, FZstd )
      Do 240 l = 1, AlnLen, 1
         FZEnt(l) = ( fe(l) - FZav ) / FZstd
  240    continue
C
      Do 1000 ng = 1, GrpKnt, 1
         If( ng .eq. NotID )    GoTo 1000
         Do 300 l = 1, AlnLen, 1
            ge(l) = Cross(l,1,ng) + Cross(l,2,ng)
  300       continue
         Call Rindex( AlnLen, ge, GRank )
         Do 320 l = 1, LNoGap, 1
            lm = AlnLen + 1 - l
            gp(l) = ge(GRank(lm))
  320       continue
         Call AveVar( gp, LNoGap, GZav, GZstd )
         Do 340 l = 1, AlnLen, 1
            GZEnt(l,ng) = ( ge(l) - GZav ) / GZstd
            If( MaxGZ(l) .lt. GZEnt(l,ng) )    MaxGZ(l) = GZEnt(l,ng)
  340       continue
C
      Write( OutFl, 5 )  GrpNam(ng), SqName( InvGrp(1,ng) ), FZav,
     *                   FZstd, GZav, GZstd
C      Write( GrupFl, 8 )
C      sn = InvGrp( 1, ng )
C      Do 360 l = 1, AlnLen, 1
C         Write( GrupFl, 9) l, Fe(l), fz(l), fp(l), Ge(l), gz(l), gp(l),
C     *                     SqPrnt( Seq(LsCum(sn)+l) ), GapMap(l)
C  360    continue
C
 1000    continue
C
      Return
      End
C
      Subroutine GeneDoc( GrupFl, MaxSqP, MxGrpP, NSP, TotSqP,
     *                    Gap, Seq, LsCum, AlnLen, LGpNam, Entropy,
     1                    Cross, GrpNam, SqPrnt, LNoGap, GrpKnt,
     2                    InvGrp, SqName, GapMap, NotID, GZEnt, FZEnt )
      Implicit None
C
C ***   Passed Variable Declarations
C
      Integer        GrupFl, MaxSqP, MxGrpP, NSP, TotSqP, Gap, AlnLen,
     *               LNoGap, GrpKnt, NotID
      Integer        Seq( TotSqP ), LsCum( NSP ), LGpNam( MxGrpP ),
     *               GapMap( MaxSqP ), InvGrp( NSP, MxGrpP )
      Real*4         Entropy( MaxSqP ), Cross(MaxSqP, 2, MxGrpP ),
     *               GZEnt(MaxSqP,MxGrpP), FZEnt(MaxSqP)
      Character*1    SqPrnt( 25 )
      Character*15   SqName( NSP )
      Character*30   GrpNam( MxGrpP )
C
C ***   Local Variable Declarations
C
C ***   Note:  The local parameter LMxSqP should be at least as large
C              as the Global parameter MaxSqP.
C
      Integer        LMxSqP
      Parameter  (   LMxSqP = 1000 )

C
C ***   Ent3V( 1=Q, 2=F, 3=f, 4=G, 5=g, 6=blank )
C ***   EntFV( 1=1, 2=2, 3=3, 4=blank )
C ***   EntGV( 1=1, 2=2, 3=3, 4=blank )
C ***   EntTV( 1=1, 2=2, 3=3, 4=blank )
C
      Integer        i, ng, sn, l, lm, Count, Ntop, NSqRes
      Integer        GRank( LMxSqP ), FRank( LMxSqP ), Ent3V( 6 ),
     *               EntFV( 4 ), EntGV( 4 ), EntTV( 4 )
      Real*4         Top10, Top20, Top30
      Real*4         gz( LMxSqP ), gp( LMxSqP )
      Character*34   EFile
      Character*60   Seq60, Ent3, EntG, EntF, TopEnt
C
    1 Format(/'Sequence|', A60 )
    2 Format( 'Entropy3|', A60 )
    3 Format( 'EntropyG|', A60 )
    4 Format( 'EntropyF|', A60 )
    5 Format( '//           Group: ', A30,
     *       /'// Master Sequence: ', A15,
     1       /'// Average Family Entropy =', F8.3, ',  Standard',
     2        ' Deviation =', F8.3,
     3       /'//  Average Group Entropy =', F8.3, ',  Standard',
     4        ' Deviation =', F8.3 )
    6 Format( 'EntropyT|', A60 )
C    7 Format(' Ntop = ', I5, ', Top10, Top20, Top30 =', 3F8.3 )
C    8 Format('   L    Fe    fz    fp      Ge    gz    gp seq  GapMap')
C    9 Format( I5, 6F8.3, 2X, A1, I8 )
   10 Format(//'//Entropy3:  Q''s=', I5, ',  F''s=', I5, ',  f''s=',I5,
     *        ',  G''s=',I5, ',  g''s=', I5 )
   11 Format('//EntropyG:  1''s=', I5, ',  2''s=', I5, ',  3''s=', I5 )
   12 Format('//EntropyF:  1''s=', I5, ',  2''s=', I5, ',  3''s=', I5 )
   13 Format('//EntropyT:  1''s=', I5, ',  2''s=', I5, ',  3''s=', I5 )
   14 Format('//There are',I6, ' residues in sequence: ', A15 )
C
C
      Do 1000 ng = 1, GrpKnt, 1
         If( ng .eq. NotID )    GoTo 1000
         Do 200 l = 1, AlnLen, 1
            gz(l) = GZEnt(l,ng)
  200       continue
         EFile = '                                  '
         EFile(1:LGpNam(ng)) = GrpNam(ng)(1:LGpNam(ng))
         EFile(LGpNam(ng)+1:LGpNam(ng)+4) = '.stu'
         Open( Unit = GrupFl, File = EFile, Status = 'New',
     *         Access = 'Sequential', Form = 'Formatted')
c     *         carriageControl = 'List' )
C
         Do 260 l = 1, 4, 1
            Ent3V(l) = 0
            EntFV(l) = 0
            EntGV(l) = 0
            EntTV(l) = 0
  260       continue
         Ent3V(5) = 0
         Ent3V(6) = 0
         NSqRes = 0
         sn = InvGrp( 1, ng )
         Call Rindex( AlnLen, gz, GRank )
         Do 320 l = 1, LNoGap, 1
            lm = AlnLen + 1 - l
            gp(l) = gz(GRank(lm))
  320       continue
         Ntop = ( LNoGap + 5 ) / 10
         Top10 = gp(NTop)
         Top20 = gp(2*NTop)
         Top30 = gp(3*NTop)
C
         Count = 0
         Seq60 = ' '
         Ent3 = ' '
         EntG = ' '
         EntF = ' '
         TopEnt = ' '
C
         Do 500 l = 1, AlnLen, 1
            If( Seq(LsCum(sn)+l) .eq. Gap )    GoTo 500
            Count = Count + 1
            NSqRes = NSqRes + 1
            Seq60(Count:Count) = SqPrnt( Seq(LsCum(sn)+l) )
            If( GapMap(l) .ge. 1 )    Then
               If( FZEnt(l) .ge. 4.0 )    Then
                  Ent3(Count:Count) = 'F'
                  Ent3V(2) = Ent3V(2) + 1
                  EntF(Count:Count) = '1'
                  EntFV(1) = Ent3V(1) + 1
               ElseIf( FZEnt(l) .ge. 3.0 )    Then
                  Ent3(Count:Count) = 'f'
                  Ent3V(3) = Ent3V(3) + 1
                  EntF(Count:Count) = '2'
                  EntFV(2) = Ent3V(2) + 1
               ElseIf( FZEnt(l) .ge. 2.0 )    Then
                  EntF(Count:Count) = '3'
                  EntFV(3) = Ent3V(3) + 1
               Else
               EndIf
               If( gz(l) .ge. 4.0 )    Then
                  Ent3(Count:Count) = 'G'
                  Ent3V(4) = Ent3V(4) + 1
                  EntG(Count:Count) = '1'
                  EntGV(1) = EntGV(1) + 1
               ElseIf( gz(l) .ge. 3.0 )    Then
                  Ent3(Count:Count) = 'g'
                  Ent3V(5) = Ent3V(5) + 1
                  EntG(Count:Count) = '2'
                  EntGV(2) = EntGV(2) + 1
               ElseIf( gz(l) .ge. 2.0 )    Then
                  EntG(Count:Count) = '3'
                  EntGV(3) = EntGV(3) + 1
               Else
               EndIf
               If( gz(l) .ge. Top10 )    Then
                  TopEnt(Count:Count) = '1'
                  EntTV(1) = EnTGV(1) + 1
               ElseIf( gz(l) .ge. Top20 )    Then
                  TopEnt(Count:Count) = '2'
                  EntTV(2) = EnTGV(2) + 1
               ElseIf( gz(l) .ge. Top30 )    Then
                  TopEnt(Count:Count) = '3'
                  EntTV(3) = EnTGV(3) + 1
               Else
               EndIf
               If( FZEnt(l) .ge. 3.0   .and.   gz(l) .ge. 3.0 )    Then
                  Ent3(Count:Count) = 'Q'
                  Ent3V(1) = Ent3V(1) + 1
                  Ent3V(2) = Ent3V(2) - 1
                  Ent3V(3) = Ent3V(3) - 1
                  Ent3V(4) = Ent3V(4) - 1
                  Ent3V(5) = Ent3V(5) - 1
               Else
               EndIf
            Else
            EndIf
            If( Count .ge. 60 )    Then
               Write( GrupFl, 1 )  Seq60
               Write( GrupFl, 2 )  Ent3
               Write( GrupFl, 3 )  EntG
               Write( GrupFl, 4 )  EntF
               Write( GrupFl, 6 )  TopEnt
               Count = 0
               Seq60 = ' '
               Ent3 = ' '
               EntG = ' '
               EntF = ' '
               TopEnt = ' '
            Else
            EndIf
  500       continue
         If( Count .gt. 0 )    Then
            Write( GrupFl, 1 )  Seq60
            Write( GrupFl, 2 )  Ent3
            Write( GrupFl, 3 )  EntG
            Write( GrupFl, 4 )  EntF
            Write( GrupFl, 6 )  TopEnt
         Else
         EndIf
         Write( GrupFl, 10 )  ( Ent3V(i), i = 1, 5 )
         Write( GrupFl, 11 )  ( EntFV(i), i = 1, 3 )
         Write( GrupFl, 12 )  ( EntGV(i), i = 1, 3 )
         Write( GrupFl, 13 )  ( EntTV(i), i = 1, 3 )
         Write( GrupFl, 14 )  NSqRes, SqName( InvGrp(1,ng) )
         Close( Unit = GrupFl, Status = 'Keep' )
 1000    continue
C
      Return
      End
C
C
      SUBROUTINE RINDEX(N,Array,INDX)
C
C ***   Modified by J. Burkardt from "Numerical Recipes" Real index routine.
C
      DIMENSION ARRAY(N),INDX(N)
C
      DO 11 J=1,N
        INDX(J)=J
11      CONTINUE
      L=N/2+1
      ir=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          q=Array(INDXT)
        ELSE
          INDXT=INDX(ir)
          q=Array(INDXT)
          INDX(ir)=INDX(1)
          ir=ir-1
          IF(ir.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.ir)THEN
          IF(J.LT.ir)THEN
            IF(Array(INDX(J)).LT.Array(INDX(J+1)))J=J+1
          ENDIF
          IF(q.LT.Array(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=ir+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
C
C      Subroutine RMotif( Termnl, Keybrd, MaxSqP, NMotif, IMotif,
C     *                   LMotif, BMotif, KMotif, KMPos )
C
C ***   Reads a file that specifies which sequence positions are part of
C ***   motifs, sequence regions that have a high concentration of highly
C ***   conserved residues that are likely to be important for maintaning
C ***   the structure and function of the macromolecular sequence family.
C ***   Motifs are specified by their index or column number in the
C ***   multiple sequence alignment being analyzed.  Each motif is
C ***   represented separately in the file.  Each motif entry begins with
C ***   a line in the file on which the word motif is the first item.
C ***   Further items in the entry can be on the same line or on
C ***   subsequent lines.  Following the entry beginning of the word motif
C ***   the alignment position numbers or indices follow as integers
C ***   separated by either spaces or commas.  The entry is terminated
C ***   by a double slash, //.  Thus a complete entry could be:
C
C    motif   23, 24, 25, 26, 27, 28  //
C
C ***   Note:  the alignment position numbers for a motif need not be
C ***          consecutive.
C
C      Implicit None
C
C ***   Passed Variable and Parameter Declarations
C
C      Integer     Termnl, KeyBrd, MaxSqP, KMotif, KMPos
C      Integer     NMotif( MaxSqP ), IMotif( MaxSqP ), LMotif( MaxSqP ),
C     *            BMotif( MaxSqP )
C
C ***   Local Variable and Parameter Declarations
C
C ***   Variables used in parsing input lines
C ***   LMxSqP should be set to at least as large as MaxSqP that is
C ***   passed in from the main program
C
C      INTEGER       NumbrP, MaxLin, LMxSqP, MtUnit
C
C      PARAMETER  (  NumbrP = 30,   MaxLin = 80, LMxSqP = 1000,
C     *              MtUnit = 3  )
C
C      INTEGER       Count, IntNum, IStyle
C      INTEGER       Digits( NumbrP ), LayOut( NumbrP )
C      CHARACTER*31  TmpStr, Numbrs( NumbrP )
C      CHARACTER*80  ILine, JLine, UCLine        !   must be character*MaxLin
C
C      Integer       i, n, j, NBlock, nb, ib, ibm, ie, MFile, KPos
C      Character*10  IntTmp
C      Character*80  MFName
C      Logical       empty
C      Save          IStyle
C      Data   IStyle  /  1  /
C
C
C    1 Format(//' Please enter the name of the file of motif',
C     *         ' descriptions' )
C    2 Format( A80 )
C    3 Format( I10 )
C    4 Format(/' Motif', I4, ',  Length=', I4, ',   begins at', I4,       DeBug
C     *        ' of array IMotif,'                                        DeBug
C     1       /' and contains the alignment positions listed below.' )    DeBug
C    5 Format('  ', 15( I4:, ',' ) )                                      DeBug
C
C      MFile = MtUnit
C      KMotif = 0
C      KMPos = 0
C      Do 150 i = 1, MaxSqP, 1
C         NMotif(i) = 0
C         IMotif(i) = 0
C         LMotif(i) = 0
C         BMotif(i) = 0
C  150    continue
C      Write( Termnl, 1 )
C      Read( KeyBrd, 2 )  MFName
C      Call PakNam( MFName, 80, empty )
C      If( empty )    Return
C      Open( Unit = MFile, File = MFName, Status = 'old' )
C
C ***   Look for the beginning of a motif entry, the word "motif".
C
C  200 Read( MFile, 2, End = 500 )  ILine
C      Call PakNam( ILine, MaxLin, empty )
C      If( empty )    GoTo 200
C      UCLine = ILine
C      Call AllCap( UCLine, MaxLin )
C      If( UCLine(1:5) .ne. 'MOTIF' )    GoTo 200
C  220 KMotif = KMotif + 1
C      BMotif(KMotif) = KMPos + 1
C      LMotif(KMotif) = 0
C  250 Call NumLst( Termnl, KeyBrd, Count, NumbrP, MaxLin, UCLine,
C     *             Numbrs, Digits, LayOut )
C      Do 300 n = 1, Count, 1
C         If( LayOut(n) .eq. IStyle )    Then
C            IntTmp = '          '
C            IntTmp(11-Digits(n):10) = Numbrs(n)(1:Digits(n))
C            Read( IntTmp, 3 )  j
C            NMotif(j) = KMotif
C            KMPos = KMPos + 1
C            IMotif(KMPos) = j
C            LMotif(KMotif) = LMotif(KMotif) + 1
C         Else
C            If( Numbrs(n)(1:2) .eq. '//' )    GoTo 200
C         EndIf
C  300    continue
C
C  400 Read( MFile, 2, End = 500 )  ILine
C      Call PakNam( ILine, MaxLin, empty )
C      If( empty )    GoTo 400
C      UCLine = ILine
C      Call AllCap( UCLine, MaxLin )
C      If( UCLine(1:5) .eq. 'MOTIF' )    Then
C         GoTo 220
C      Else
C         GoTo 250
C      EndIf
C
C  500 Close( Unit = MFile, status = 'keep' )
C
C      Do 700 n = 1, KMotif, 1                                            DeBug
C         Write( Termnl, 4 )  n, LMotif(n), BMotif(n)                     DeBug
C         NBlock = ( LMotif(n) + 14 ) / 15                                DeBug
C         ib = BMotif(n) - 15                                             DeBug
C         Do 650 nb = 1, NBlock, 1                                        DeBug
C            ib = ib + 15                                                 DeBug
C            ibm = ib - 1                                                 DeBug
C            ie = ib + 14                                                 DeBug
C            If( ie .gt. ibm + LMotif(n) )    ie = ibm + LMotif(n)        DeBug
C            Write( Termnl, 5 )  ( IMotif(i), i = ib, ie, 1 )             DeBug
C  650       continue                                                     DeBug
C  700    continue                                                        DeBug
C
C      Return
C      End
C
C      Subroutine WrVol( MaxSqP, TotSqP, NSP, NS, MxSymP, AlnLen, Gap,
C     *                  LNoGap, LsCum, Seq, GrupId, NoGaps, GapMap,
C     1                  NSqUsed, FrGap, BackFr, ASize, UnKnown,
C     2                  NMotif, Entropy, Consen, BestSqT, sqprnt,
C     3                  MaxVol, MinVol, AveVol, SigVol )
C
C      Implicit None
C
C ***   Passed variable declarations
C
C      Integer      MaxSqP, TotSqP, NSP, NS, MxSymP, AlnLen, Gap,
C     *             LNoGap, NSqUsed, OutFl, ASize, UnKnown
C      Integer      LsCum( NSP ), Seq( TotSqP ), NoGaps( MaxSqP ),
C     *             GapMap( MaxSqP ), GComp( 0:MxSymP ), GrupId( NSP ),
C     *             NMotif( MaxSqP )
C      Real*4       FrGap( MaxSqP ), BackFr( MxSymP ), MaxVol( MaxSqP ),
C     *             MinVol( MaxSqP ), AveVol( MaxSqP ), SigVol( MaxSqP),
C     *             Entropy(MaxSqP)
C      Character*1  BestSqT(MaxSqP), Consen(MaxSqP), SqPrnt(25)
C
C ***   Local Variable declarations
C
C      Integer     VlUnit, GpUnit, SqUnit
C      Parameter  ( VlUnit = 11, GpUnit = 12, SqUnit = 13 )
C
C      Integer     VolFl, GapFl, SeqFl, l, n, aa, aaknt, maxaa,
C     *            tc(0:25)
C      Real*4      fraa, frvol
C
C    1 Format(' ', I4, 2X, A1, ' ', F6.2, 2X, A1, F8.2, 3F7.1, F9.2,
C     *       F6.2, I5, I4 )
C    2 Format('//   Position on row     --     Item description.',
C     *      /'//'
C     1      /'//   1   --   Alignment Index (column in alignment)',
C     2      /'//   2   --   Consensus, most common, residue symbol',
C     3      /'//   3   --   Fraction of consensus residue in column',
C     4      /'//   4   --   Highest entropy residue',
C     5      /'//   5   --   Entropy for alignment column',
C     6      /'//   6   --   Average volume of residues in column',
C     7      /'//   7   --   Maximum volume of residues in column',
C     8      /'//   8   --   Minimum volume of residues in column',
C     9      /'//   9   --   Standard deviation of the volume of',
C     A                      ' residues in column',
C     B      /'//  10   --   ( Max Vol - Min Vol ) / Ave Vol',
C     C      /'//  11   --   Number of gaps (indels) in column',
C     D      /'//  12   --   Motif identification number',
C     E      /'//' )
C    3 Format(' ', I4, 2X, A1, ' ', F6.2, 2X, 3F7.1, F9.2, F6.2, I5 )
C    4 Format('//   Position on row     --     Item description.',
C     *      /'//'
C     1      /'//   1   --   Alignment Index (column in alignment)',
C     2      /'//   2   --   Consensus, most common, residue symbol',
C     3      /'//   3   --   Fraction of consensus residue in column',
C     4      /'//   4   --   Average volume of residues in column',
C     5      /'//   5   --   Maximum volume of residues in column',
C     6      /'//   6   --   Minimum volume of residues in column',
C     7      /'//   7   --   Standard deviation of the volume of',
C     8                      ' residues in column',
C     9      /'//   8   --   ( Max Vol - Min Vol ) / Ave Vol',
C     A      /'//   9   --   Number of gaps (indels) in column',
C     B      /'//' )
C    5 Format(' ', I4, 2X, A1, F6.2, 20( I4, A1: ) )
C    6 Format('//   Position on row     --     Item description.',
C     *      /'//'
C     1      /'//   1   --   Alignment Index (column in alignment)',
C     2      /'//   2   --   Consensus, most common, residue symbol',
C     3      /'//   3   --   Fraction of consensus residue in column',
C     4      /'//   4->      Symbol and number of occurrences in column',
C     5      /'//' )
C
C      VolFl = VlUnit
C      GapFl = GpUnit
C      SeqFl = SqUnit
C      Open( Unit = VolFl, File = 'NoGaps.vol', Status = 'new',
C     *      Access = 'Sequential', Form = 'formatted',
C     *      CarriageControl = 'List' )
C      Open( Unit = GapFl, File = 'GapPos.vol', Status = 'new',
C     *      Access = 'Sequential', Form = 'formatted',
C     *      CarriageControl = 'List' )
C      Open( Unit = SeqFl, File = 'Conserved.pos', Status = 'new',
C     *      Access = 'Sequential', Form = 'formatted',
C     *      CarriageControl = 'List' )
C      Write( VolFl, 2 )
C      Write( GapFl, 4 )
C      Write( SeqFl, 6 )
C
C      Do 400 l = 1, AlnLen, 1
C         Do 200 n = 1, 25, 1
C            tc(n) = 0
C  200       continue
C         Do 300 n = 1, NS, 1
C            If( GrupId(n) .ge. 1 )    Then
C               aa = Seq( LsCum(n) + l )
C               tc( aa ) = tc( aa ) + 1
C            Else
C            EndIf
C  300       continue
C         maxaa = 0
C         aaknt = 0
C         Do 320 aa = 1, ASize, 1
C            aaknt = aaknt + tc(aa)
C            If( tc(aa) .gt. maxaa )    maxaa = tc(aa)
C  320       continue
C         fraa = Float( maxaa ) / Float( aaknt )
C         frvol = ( MaxVol(l) - MinVol(l) ) / AveVol(l)
C         If( tc(gap) .eq. 0 )    Then
C            Write( VolFl, 1 )  l, Consen(l), fraa, BestSqT(l),
C     *                         Entropy(l), AveVol(l), MaxVol(l),
C     1                         MinVol(l), SigVol(l), FrVol, tc(gap),
C     2                         NMotif(l)
C         Else
C            If( aaknt .gt. 3 )    Write( GapFl, 3 )  l, Consen(l),
C     1                                 fraa, AveVol(l), MaxVol(l),
C     2                                 MinVol(l), SigVol(l), FrVol,
C     3                                 tc(gap)
C         EndIf
C         Write( SeqFl, 5 )  l, Consen(l), fraa, ( tc(aa), sqprnt(aa),
C     *                                            aa = 1, ASize, 1 )
C  400    continue
C
C      Close( Unit = VolFl, Status = 'keep' )
C      Close( Unit = GapFl, Status = 'keep' )
C      Close( Unit = SeqFl, Status = 'keep' )
C      Return
C      End
C
      Subroutine WrCons( SeqFl, AlnLen, GrpKnt, PFLen, PFile, MaxSqP,
     *                   MxGrpP, Consen, GCons, BestSqT, BestSq,
     1                   GrpNam )
      Implicit None
C
C ***   Passed variable declararations
C
      Integer        SeqFl, AlnLen, GrpKnt, PFLen, MaxSqP, MxGrpP
      Character*1    Consen( MaxSqP ), GCons( MaxSqP, 2, MxGrpP ),
     *               BestSqT( MaxSqP ), BestSq( MaxSqP, 2, MxGrpP )
      Character*24   PFile
      Character*30   GrpNam
C
C ***   Local varaible declarations
C
      Integer        i, g, ib, ie, nb, Blocks
C
    1 Format( '>', A30 )
    2 Format( '>', A24 )
    3 Format( ' ', 7( ' ', 10(A1:) ))
C
      Blocks = ( AlnLen + 69 ) / 70
      Open( Unit = SeqFl, File = 'mostcommon.aa', Status = 'new',
     *      Form = 'formatted',
     *      Access = 'sequential' )
C     *      CarriageControl = 'list', Form = 'formatted',
      ib = -69
      Write( SeqFl, 2 )  PFile
      Do 100 nb = 1, Blocks, 1
         ib = ib + 70
         ie = ib + 69
         If( ie .gt. AlnLen )    ie = AlnLen
         Write( SeqFl, 3 )  ( Consen(i), i = ib, ie, 1 )
  100    continue
      Do 200 g = 1, GrpKnt, 1
         ib = -69
         Write( SeqFl, 1 )  GrpNam( g )
         Do 150 nb = 1, Blocks, 1
            ib = ib + 70
            ie = ib + 69
            If( ie .gt. AlnLen )    ie = AlnLen
            Write( SeqFl, 3 )  ( GCons(i,1,g), i = ib, ie, 1 )
  150       continue
  200    continue
      Close( Unit = SeqFl, Status = 'keep' )
C
      Open( Unit = SeqFl, File = 'highentropy.aa', Status = 'new',
     *      Form = 'formatted',
     *      Access = 'sequential' )
c     *      CarriageControl = 'list', Form = 'formatted',
      ib = -69
      Write( SeqFl, 2 )  PFile
      Do 300 nb = 1, Blocks, 1
         ib = ib + 70
         ie = ib + 69
         If( ie .gt. AlnLen )    ie = AlnLen
         Write( SeqFl, 3 )  ( BestSqT(i), i = ib, ie, 1 )
  300    continue
      Do 400 g = 1, GrpKnt, 1
         ib = -69
         Write( SeqFl, 1 )  GrpNam( g )
         Do 350 nb = 1, Blocks, 1
            ib = ib + 70
            ie = ib + 69
            If( ie .gt. AlnLen )    ie = AlnLen
            Write( SeqFl, 3 )  ( BestSq(i,1,g), i = ib, ie, 1 )
  350       continue
  400    continue
      Close( Unit = SeqFl, Status = 'keep' )
      Return
      End
