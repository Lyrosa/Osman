import time
import re
import mysql.connector
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def main():
#==============================================================================
    file1 = open("blastresults.txt", "w") .close()
#==============================================================================
#   file1= open ("blastresulttest.txt", "w")
#==============================================================================
    file1 = open("blastresults.txt", "a")
#==============================================================================
     
     
    bestand = open("seqchamp.txt", "r")
    leeg = ""
    for regel in bestand:
        leeg += regel
    s= leeg.split("@HWI")
    s.remove("")
    
    d = maakDict(s)
    
    count = 1
    
    for dkey, dvalue in d.items():
        

        if count < 5:

            nrblast= ("\n\n\n\n-----------------blast number: "+ str(count)+"-----------------------\n\n\n")
            head = ("@HWI"+dkey)
            head = ("header: "+head)
            headseq = ("sequence: "+ dvalue)
            
            print(nrblast)
            print(head)
            print(headseq)
    
            file1.write ("\n"+str(nrblast))
            file1.write ("\n"+str(head))
            file1.write ("\n"+str(headseq))
    
            resultaat = blast(dvalue, file1, dkey)
            mainUpdate(resultaat, dkey, dvalue)
            if resultaat == False:
                count += 1
                print("No results")
            else: 
                endresult = ("----------------------end--------------------")
                print (endresult)
                count +=1
#==============================================================================
#     file1.close()
#==============================================================================
    
def maakDict(s):
    dictionary= {}
    for i in s:
       nlijst = i.split("\t")
       h = nlijst[0]
       seq = nlijst[1]
       dictionary.update({h:seq})
    return (dictionary)
    
    
def mainUpdate(bool, head, seq):
    conn = mysql.connector.connect(host="ithurtswhenip.nl",user = "pg3", password = "pg3", db = "pg3", port = 3307)
    cursor = conn.cursor()

    query1 = ("INSERT INTO `pg3`.`main_seq`(`HWI_code`, `seq`, `match_found`) VALUES ('"+str(head)+"', '"+str(seq)+"', "+str(bool)+");")

#==============================================================================
#     query2 = ("INSERT INTO `pg3`.`results_blasts`(`result1_accession_code_result`, `result1_length_result`, `result1_identity`, `result1_positives`, `result1_E_value`, `result1_gaps`, `result1_query`, `result1_match`, `result1_subject`) VALUES ('"+str(accession)+"',"+str(length)+",'"+str(identity)+"','"+str(positives)+"',"+str(e_value)+","+str(gaps)+", '"+str(query)+"','"+str(match)+"''"+str(subject)+"');")
#==============================================================================
    cursor.execute(query1)
    print (query1)
    conn.commit()
    cursor.close()
    conn.close()  

        



def blast(sequence, file1, head):
    #try:
    bool = True
    count = 1
    results_handle = NCBIWWW.qblast("blastx","nr", sequence, hitlist_size=5, expect=1)
#==============================================================================
#     conn = mysql.connector.connect(host = "ithurtswhenip.nl",user = "richard", password = "richard", db = "pg3", port = 3307) 
#     cursor = conn.cursor ()
#==============================================================================
##        bestand= open ("blast_report.xml", "w")
##        bestand.writelines (results_handle.readlines())
##        bestand.close()
##        result = open("blast_report.xml", "r")
        #blast_records = NCBIXML.parse(result)
        #blast_record = next(blast_records)
##        
##    except TimeoutError:
##        print ("TimeoutError: Hervat blasten over 60 seconden.")
##        for i in range (60):
##            print (60-i)
##            time.sleep(1)
##            break
##    except ConnectionAbortedError:
##        print ("ConnectionAbortedError: Poging tot hervatten van blasten over 60 seconden.")
##        for i in range (60):
##            print (60-i)
##            time.sleep(1)
##            break
##    except:
##        print("Onbekende Error: Poging tot hervatten in 60 seconden.")
##        for i in range (60):
##            print (60-i)
##            time.sleep(1)
##            break
#==============================================================================
#         result = open("blast_report.xml", "r")
#         blast_records = NCBIXML.parse(result)
#         blast_record = next(blast_records)
#==============================================================================
    
    #bijgekomen is Accesion, eiwit score identity en positives
    blast_results = results_handle
    blast_records = NCBIXML.parse(results_handle)

    E_VALUE_THRESH = 1
    for blast_record in blast_records:
        desc = blast_record.descriptions
        if desc == []:
            bool = False
        else:
            numhits = len(desc)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    nralignment = ('\n------Alignment'+str(count)+'-------\n')
                    header = alignment.title
                    length = alignment.length
                    accession = alignment.accession
                    eiwit  =   alignment.hit_def
                    score  = hsp.score
                    e_value = hsp.expect
                    identity = hsp.identities
                    positives = hsp.positives
                    gaps    =  hsp.gaps
                    query   =  hsp.query
                    match   =  hsp.match
                    subject =  hsp.sbjct
                    count += 1
                    
                    
#==============================================================================
#                     query2 = ("INSERT INTO `pg3`.`results_blasts`(`result1_accession_code_result`, `result1_length_result`, `result1_identity`, `result1_positives`, `result1_E_value`, `result1_gaps`, `result1_query`, `result1_match`, `result1_subject`) VALUES ('"+str(accession)+"',"+str(length)+",'"+str(identity)+"','"+str(positives)+"',"+str(e_value)+","+str(gaps)+", '"+str(query)+"','"+str(match)+"''"+str(subject)+"');")
#                     query3 = ("INSERT INTO `pg3`.`results_blasts`(`result2_accession_code_result`, `result2_length_result`, `result2_identity`, `result2_positives`, `result2_E_value`, `result2_gaps`, `result2_query`, `result2_match`, `result2_subject`) VALUES ('"+str(accession)+"',"+str(length)+",'"+str(identity)+"','"+str(positives)+"',"+str(e_value)+","+str(gaps)+", '"+str(query)+"','"+str(match)+"''"+str(subject)+"');")
#                     query4 = ("INSERT INTO `pg3`.`results_blasts`(`result3_accession_code_result`, `result3_length_result`, `result3_identity`, `result3_positives`, `result3_E_value`, `result3_gaps`, `result3_query`, `result3_match`, `result3_subject`) VALUES ('"+str(accession)+"',"+str(length)+",'"+str(identity)+"','"+str(positives)+"',"+str(e_value)+","+str(gaps)+", '"+str(query)+"','"+str(match)+"''"+str(subject)+"');")
#                     query5 = ("INSERT INTO `pg3`.`results_blasts`(`result4_accession_code_result`, `result4_length_result`, `result4_identity`, `result4_positives`, `result4_E_value`, `result4_gaps`, `result4_query`, `result4_match`, `result4_subject`) VALUES ('"+str(accession)+"',"+str(length)+",'"+str(identity)+"','"+str(positives)+"',"+str(e_value)+","+str(gaps)+", '"+str(query)+"','"+str(match)+"''"+str(subject)+"');")
#                     query6 = ("INSERT INTO `pg3`.`results_blasts`(`result5_accession_code_result`, `result5_length_result`, `result5_identity`, `result5_positives`, `result5_E_value`, `result5_gaps`, `result5_query`, `result5_match`, `result5_subject`) VALUES ('"+str(accession)+"',"+str(length)+",'"+str(identity)+"','"+str(positives)+"',"+str(e_value)+","+str(gaps)+", '"+str(query)+"','"+str(match)+"''"+str(subject)+"');")
#==============================================================================

                    if count == 1:
                        conn = mysql.connector.connect(host="ithurtswhenip.nl",user = "pg3", password = "pg3", db = "pg3", port = 3307)
                        cursor = conn.cursor()
                        query2 = ("INSERT INTO `pg3`.`results_blasts`(`result1_accession_code_result`, `result1_length_result`, `result1_identity`, `result1_positives`, `result1_E_value`, `result1_gaps`, `result1_query`, `result1_match`, `result1_subject`, `main_seq_HWI_code`) VALUES ('"+str(accession)+"',"+str(length)+",'"+str(identity)+"','"+str(positives)+"',"+str(e_value)+","+str(gaps)+", '"+str(query)+"','"+str(match)+"','"+str(subject)+"', '"+str(header)+"');")
                        cursor.execute(query2)
                        print(query2)
                        conn.commit()
                        cursor.close()
                        conn.close() 
                        count += 1
                    elif count < 6:
                        conn = mysql.connector.connect(host="ithurtswhenip.nl",user = "pg3", password = "pg3", db = "pg3", port = 3307)
                        cursor = conn.cursor()
                        query3 = ("UPDATE `pg3`.`results_blast` SET ` result"+str(count)+"_accession_code_result`='"+str(accession)+"', `result"+str(count)+"_length_result`="+str(length)+", `result"+str(count)+"_identity`='"+str(identity)+"', `result"+str(count)+"_positives`='"+str(positives)+"', `result"+str(count)+"_E_value`="+str(e_value)+", `result"+str(count)+"_gaps`="+str(gaps)+", `result"+str(count)+"_query`='"+str(query)+"', `result"+str(count)+"_match`='"+str(match)+"', `result"+str(count)+"_subject`='"+str(subject)+"' WHERE `main_seq_HWI_code`='"+str(header)+"';")
                        cursor.execute(query3)
                        print(query3)
                        conn.commit()
                        cursor.close()
                        conn.close()
                        count += 1
                                               



                    
                    print (nralignment)
                    print ('sequence :', header)
                    print ('length: ', length)
                    print ('Accession: ',accession )
                    print ('Eiwit: ', eiwit)
                    print ('Score: ',score)
                    print ('e value:',e_value)
                    print ('Identity: ', identity)
                    print ('Positives: ',positives)
                    print ('Gaps:   ',gaps)
                    print (query)
                    print (match)
                    print (subject)
                    
                    file1.write ("\n"+str(nralignment))
                    file1.write ("\n Sequence titel:  "+str(header))
                    file1.write ("\n Length:    "+str(length))
                    file1.write ("\n Accession:  "+str(accession))
                    file1.write ("\n Eiwit:      "+str(eiwit))
                    file1.write ("\n Score:      "+str(score))
                    file1.write ("\n E_value:    "+str(e_value))
                    file1.write ("\n Identity%:   "+str(identity))
                    file1.write ("\n Positives:  "+str(positives))
                    file1.write ("\n Gaps:       "+str(gaps))
                    file1.write ("\n"+str(query))
                    file1.write ("\n"+str(match))
                    file1.write ("\n"+str(subject))
    return (bool)            
            
                
                

main()
