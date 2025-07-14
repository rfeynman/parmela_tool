'''parmela tool, used to optimize the beamline.
scan the value to get min emittance at different bunch length
'''
import os
import pandas as pd
import numpy as np

parmela='C:/LANL/parmela.exe '
# parmela = 'wine ~/.wine/drive_c/LANL/parmela.exe '


def getvar(filename):
    file = open(filename, 'rU')
    lines = file.readlines()
    file.close()
    mark = []
    tempmark = []
    step = []
    tempstep = []
    left_range = []
    templeft_range = []
    right_range = []
    tempright_range = []
    pos = []
    for line in lines:
        line = line.strip()
        var = line.split()
        ln = len(var)
# when open file with 'rU' doesn't work (ASCII decode doesn't work)
# change 'rU' to 'rb' and use following loop to see which word can't decode
#        for i in range(ln):
#            print(var[i])
#            var[i] = var[i].decode('ASCII')
        if var != []:
            if (var[0] == '!@var' and var[1] == '1'):
                tempmark.append(var[2])
                tempstep.append(var[3])
                templeft_range.append(var[4])
                tempright_range.append(var[5])
            if (var[0] == '!@subs'):
                if(ln >= 4 and var[3] in tempmark):
                    mark.append(var[3])
                    if(ln >= 6 and var[-2] == 'element'):
                        pos.append(int(var[-1]))
                if(ln >= 6 and var[5] in tempmark):
                    mark.append(var[5])
                    if(ln >= 8 and var[-2] == 'element'):
                        pos.append(int(var[-1]))
    ln1 = len(mark)
    ln2 = len(tempmark)
    for i in range(ln1):
        for j in range(ln2):
            if tempmark[j] == mark[i]:
                step.append(tempstep[j])
                left_range.append(templeft_range[j])
                right_range.append(tempright_range[j])
    return mark, step, left_range, right_range, pos


def rewriteFile(filename, mark, value):
    file = open(filename, 'rU')
    lines = file.readlines()
    file.close()
    k = 0
    for line in lines:
        line = line.strip()
        subs = line.split()
        if subs != []:
            if (subs[0] == '!@subs'):
                ln = len(subs)
                if(ln >= 4 and subs[3] == mark):
                    for i in range(int(subs[1])):
                        new_line = lines[k + i + 1]
                        new_line = new_line.split()
                        new_line[int(subs[2])] = value
                        new_line = ' '.join(new_line) + '\n'
                        lines[k + i + 1] = new_line
                if(ln >= 6 and subs[5] == mark):
                    for i in range(int(subs[1])):
                        new_line = lines[k + i + 1]
                        new_line = new_line.split()
                        new_line[int(subs[4])] = value
                        new_line = ' '.join(new_line) + '\n'
                        lines[k + i + 1] = new_line
                if(ln >= 4 and subs[3] == ('-' + mark)):
                    new_line = lines[k + i + 1]
                    new_line = new_line.split()
                    value = str(-1 * float(value))
                    new_line[int(subs[2])] = value
                    new_line = ' '.join(new_line) + '\n'
                    lines[k + i + 1] = new_line
                if(ln >= 6 and subs[5] == ('-' + mark)):
                    new_line = lines[k + i + 1]
                    new_line = new_line.split()
                    value = str(-1 * float(value))
                    new_line[int(subs[4])] = value
                    new_line = ' '.join(new_line) + '\n'
                    lines[k + i + 1] = new_line
        k += 1
    file = open(filename, 'w')
    file.writelines(lines)
    file.close()


def get_min_emittance():
    with open('EMITTANCE.TBL') as raw:
        rawdata = [line.split() for line in raw]
    raw.close()
    data = pd.DataFrame(rawdata[400:], columns=rawdata[84][1:])
    data = data.apply(pd.to_numeric)
    emit = data['Xn(mm-mrad)'].values
    min_emit = min(emit)
    return min_emit


def get_beam_size():
    with open('EMITTANCE.TBL') as raw:
        rawdata = [line.split() for line in raw]
    raw.close()
    data = pd.DataFrame(rawdata[400:], columns=rawdata[84][1:])
    data = data.apply(pd.to_numeric)
    emit = data['Xn(mm-mrad)'].values
    pos = emit.index(min(emit)) + 300

    data = pd.DataFrame(rawdata[100:pos], columns=rawdata[84][1:])
    data = data.apply(pd.to_numeric)
    size = data['Xrms(mm)'].values
    return max(size), size[-1]


def main():
    inputfilename = 'sp2.acc'
    outfilename = 'OUTPAR.TXT'
    foldername = 'scan_results'
    wfilename = 'resultfile.txt'
    mark, step, left_range, right_range, pos = getvar(inputfilename)
    emittance = []
    if mark == []:
        print('No change of the inputfile, please check the parameter')
    else:
        init_N = 1
        N = 50
        Ni = (float(right_range[2]) - float(left_range[2])) / float(step[2])
        Ni = int(Ni)
        for j in range(N - init_N + 1):
            value0 = str(float('{0:.8f}'.format((j + init_N) * float(step[0]))))
            rewriteFile(inputfilename, '0', value0)  # write the bunch length (cm)
            value1 = str(float('{0:.8f}'.format((j + init_N) * float(step[1]))))
            rewriteFile(inputfilename, '1', value1)  # write the bunch length (degree)
            bunch_length = 10 * (j + init_N)
            print('length cycle:', j,' bunch length (ps):', bunch_length)
            min_emit = []
            for i in range(Ni + 1):
                value = str(float('{0:.6f}'.format(float(left_range[2]) + i * float(step[2]))))
                rewriteFile(inputfilename, '3', value)  # write the solenoid field
                os.system(parmela + inputfilename)
                emit = get_min_emittance()
                min_emit.append(emit)
                print('field cycle:', i,' field strength (Gauss):', value, 'emittance:', emit)
            min_field = min_emit.index(min(min_emit)) * float(step[2]) + float(left_range[2])
            emittance.append([bunch_length, min_field, min(min_emit)])
            print('min emittance: (mm-mrad):', min(min_emit))
            rewriteFile(inputfilename, '3', str(min_field))
            # os.system('cp ' + inputfilename + ' ' + str(bunch_length) + '_sp.acc')
            os.system('mv EMITTANCE.TBL EMITTANCE_' + str(bunch_length) + '.TBL')
        os.system('mv EMITTANCE*.TBL ' + foldername)
        # os.system('mv *_sp.acc ' + foldername)
        emittance = np.array(emittance)
        np.savetxt('emittance.csv', emittance, delimiter=',', fmt = '%1.5f')
        os.system('mv emittance.csv ' + foldername)
        print('done')


if __name__ == '__main__':
    main()
