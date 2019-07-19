from django.shortcuts import render
from .rna_takefuji import takefuji
from .rna_olke import olke
from .rna_uv import u_v
from pathlib import Path
import os, time, re
from RNA.settings import MEDIA_ROOT
from django.core.mail import send_mail
from django.conf import settings
from django.core.files.storage import FileSystemStorage


def send_email(email, dot_bracket, seq):
    subject = "The secondary structure hence predicted"
    mes = 'The sequence you entered was: '+seq+'\n'
    for i in dot_bracket:
        mes += i
        mes += "\nLink to visualization: "+"http://nibiru.tbi.univie.ac.at/forna/forna.html?id=fasta&file=>header\\n"+seq+"\\n"+str(i)

    message = "The dot-bracket notation(s):\n" + mes +"\n\n\nThank you for letting us serve you!\n" 
    email_from = settings.EMAIL_HOST_USER
    recipient_list = [email]
    send_mail( subject, message, email_from, recipient_list)


def predict(request):
    if request.method == 'POST':
        seq = request.POST['input']
        seq = ''.join(re.split('\t| ', seq))
        keep_seq = request.POST['keep_seq']
        keep_seq = ''.join(re.split('\t| ', keep_seq))
        email = request.POST['mail']
        f1 = request.FILES.get('f1', False)
        f2 = request.FILES.get('keep_seq_file', False)
        fs = FileSystemStorage()
        if f1 is not False:
            f1_name = f1.name
            f1 = fs.save(f1_name, f1)
            path = os.path.join(MEDIA_ROOT, f1_name)
            seq = open(path, 'r').read().split('\n')[0]
            seq = ''.join(re.split('\t| ', seq))
            os.remove(path)
        
        if f2 is not False:
            f2_name = f2.name
            f2 = fs.save(f2_name, f2)
            path = os.path.join(MEDIA_ROOT, f2_name)
            keep_seq = open(path, 'r').read().split('\n')[0]
            keep_seq = ''.join(re.split('\t| ', keep_seq))
            os.remove(path)
        # home = str(Path.home())
            
        # print(len(seq), len(keep_seq))
        if len(keep_seq) == 0:
            keep_seq = "." * len(seq)
            
        if len(seq) != len(keep_seq):
            return render(request, 'app/rna.html', {'message': "Length of base sequence and .x sequence should be equal!"})
        
        bases = request.POST.getlist('base')
        method = request.POST.getlist('algo')[0]
        
        if method == 'self':
            dot_bracket, output_graphs = takefuji(seq, bases, keep_seq)
        elif method == 'olke':
            dot_bracket, output_graphs = olke(seq, bases, keep_seq)
        else:
            dot_bracket, output_graphs = u_v(seq, bases, keep_seq)

        if email != '':
            send_email(email, dot_bracket, seq)

        return render(request, 'app/rna.html', {'res': dot_bracket,
                                                'message': "",
                                                'sequence': seq.upper().replace('T','U'),
                        'dot_bracket':[(i,j) for i,j in enumerate(zip(output_graphs, dot_bracket))]})

    return render(request, 'app/rna.html')

        
