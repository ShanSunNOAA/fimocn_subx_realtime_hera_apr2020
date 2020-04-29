PRO read_tcvitals_data,infile,stormNos,stormNames

; ----- read forecast track data from file

print,'beginning subroutine read_tcvitals_data'
print, 'infile: ',infile
i=0
maxline=100
array = strarr(maxline)
storms = strarr(maxline)
stormNos = strarr(maxline)
stormNames = strarr(maxline)
openr,lun, infile,/get_lun
while not eof (lun) do begin
    line=''
    readf,lun,line
    array[i] = line
    print, 'line: ',line
    result = strsplit(line,' ',/extract)
    storms[i] = strcompress(string(result[2])+'_'+string(result[1]),/remove_all)
    stormNos[i] = strcompress(string(result[1]),/remove_all)
    stormNames[i] = strcompress(string(result[2]),/remove_all)
    print, 'in read_tcvitals storm: ',storms[i]
    print, 'in read_tcvitals stormNos: ',stormNos[i]
    print, 'in read_tcvitals stormNames: ',stormNames[i]
    i=i+1
endwhile
i=i-1
storms = storms[0:i]
stormNos = stormNos[0:i]
stormNames = stormNames[0:i]
print,'storms: ',storms
print,'in read_tcvitals_data stormNos: ',stormNos
print,'in read_tcvitals_data stormNames: ',stormNames
print, 'num elements: ',n_elements(size(storms))

free_lun,lun
return

end
