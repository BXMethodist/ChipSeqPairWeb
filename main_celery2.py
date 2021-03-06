from flask import Flask, render_template, request, session
from queryUtils import GEO_query
from search import Search
from match import Match
from setup import get_settings
from flask_celery import make_celery
from flask_mail import Mail, Message
import os

app = Flask(__name__)
app.config['CELERY_BROKER_URL'] = 'amqp://localhost//'

app.config.update(
    DEBUG = False,
    MAIL_SERVER = 'smtp.tmh.tmhs',
#    MAIL_ASCII_ATTACHMENTS = TRUE,
    MAIL_DEFAULT_SENDER= 'chenlabhoustonmethodist@gmail.com',
    MAIL_USERNAME='chenlabhoustonmethodist@gmail.com',
    MAIL_PASSWORD='R10-414D',
)

app.secret_key = 'R10-414D'

celery = make_celery(app)
mail = Mail(app)

@app.route("/")
def main():
    return render_template('index.html')

@app.route('/search')
def showSearch():
    return render_template('search.html')

@app.route('/match')
def showMatch():
    return render_template('match.html')

@app.route('/query')
def showQuery():
    return render_template('query.html')

@app.route('/display')
def showDisplay():
    return render_template('display.html')

@app.route('/double')
def showDouble():
    return render_template('double.html')

@app.route('/single')
def showSingle():
    return render_template('single.html')

@app.route('/double', methods=['POST'])
def reloadDouble():
    chr = request.form['chr']
    start = request.form['start']
    end = request.form['end']
    session['src1'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position='+chr+'%3A'+start+'-' + end
    session['src2'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position='+chr+'%3A'+start+'-' + end
    return render_template('reloadDouble.html')

@app.route('/display', methods=['POST'])
def selectFeatures():
    if 'feature1' in request.form or 'feature2' in request.form:
        _feature1 = request.form['feature1']
        _feature2 = request.form['feature2']
        if _feature1 == "" and _feature2 == "":
            session['feature1'] = 'H3K4me3'
            return render_template('single.html')
        elif _feature2 == "" or _feature1 == "":
            if _feature1 == "":
                session['feature2'] = ""
                session['feature1'] = _feature2
            else:
                session['feature2'] = ""
                session['feature1'] = _feature1
            return render_template('single.html')
        else:
            if _feature1 != _feature2:
                session['feature2'] = _feature2
                session['feature1'] = _feature1
                return render_template('double.html')
            else:
                session['feature1'] = _feature1
                session['feature2'] = ""
                return render_template('single.html')
    else:
        return render_template('display.html')

@app.route('/search',methods=['POST'])
def search():
    _searchterms = request.form['searchterms']
    _Species = request.form['Species']
    _inputEmail = request.form['inputEmail']
    species = _Species if _Species != '' else 'Homo sapiens'

    settings = get_settings()

    keywords = _searchterms.split(",")
    output_prefix = keywords[0]
    output_path = './tmp/'

    cwd = settings['Chipseq']

    CallSearch.delay(output_prefix, output_path, keywords, species, cwd, _inputEmail)

    return 'We are processing your request, results will be sent to your email'

@app.route('/match',methods=['POST'])
def match():
    _feature1 = request.form['feature1']
    _feature2 = request.form['feature2']
    _Species = request.form['Species']

    _inputEmail = request.form['inputEmail']

    settings = get_settings()

    keywords1 = _feature1.split(",")
    keywords2 = _feature2.split(",")
    output_prefix1 = keywords1[0]
    output_prefix2 = keywords2[0]
    output_path = './tmp/'

    species = _Species if _Species != '' else 'Homo sapiens'
    cwd = settings['Chipseq']

    CallMatch.delay(output_prefix1, output_prefix2, output_path, keywords1, keywords2, species, cwd, _inputEmail)

    return 'We are processing your request, results will be sent to your email'


@app.route('/query', methods=['GET', 'POST'])
def query():
    f = request.files['IDlist']
    info = [str(line) for line in f.readlines()]
    _inputEmail = request.form['inputEmail']
    id_list = []
    for line in info:
        id_list.append(line.strip())
    f.close()
    output_path = './tmp/'
    settings = get_settings()
    GSMGSE_pkl = settings['GSMGSE_pkl_path']
    GSM_SRR_pkl = settings['GSMtoSRRpkl']

    CallQuery.delay(id_list, output_path, GSMGSE_pkl, GSM_SRR_pkl, _inputEmail)
    return 'We are processing your request, results will be sent to your email'

@celery.task(name='main_celery.search')
def CallSearch(output_prefix, output_path, keywords, _Species, cwd, _inputEmail):
    result_df, samples = Search(output_prefix, output_path, keywords, _Species, cwd, _inputEmail)
    species = _Species.replace(" ", '')
    output_name = output_prefix+"_"+ species + '.csv'
    send_email(_inputEmail, [result_df], output_path, [output_name])
    return

@celery.task(name='main_celery.match')
def CallMatch(output_prefix1, output_prefix2, output_path, keywords1, keywords2, species, cwd, _inputEmail):
    result_dfs = Match(output_prefix1, output_prefix2, output_path, keywords1, keywords2, species, cwd, _inputEmail)
    output_names = [output_prefix1+'.csv', output_prefix2+'.csv', output_prefix1 + "_" + output_prefix2 +'.csv']
    send_email(_inputEmail, result_dfs, output_path, output_names)
    return

@celery.task(name='main_celery.query')
def CallQuery(id_list, output_path, GSMGSE_pkl, GSM_SRR_pkl, email):
    result_df = GEO_query(id_list, output_path, GSMGSE_pkl, GSM_SRR_pkl, email)
    output_name = 'query.txt'
    send_email(email, [result_df], output_path, [output_name])
    return

def send_email(email, results, output_path, output_names):
    subject = 'ChIPSeqPair Result'
    message = 'Hello, <br> <br>' \
              'Here is your result from ChipSeqPair <br> <br>' \
              'Thanks!  <br> <br>' \
              'Chen Lab'
    msg = Message(subject=subject, html=message, recipients=[email])
    for i in range(len(results)):
        table = results[i]
        output_name = output_names[i]
        table.to_csv(output_path+output_name,encoding='utf-8')
        with app.open_resource(output_path+output_name) as table_file:
            msg.attach(output_name, output_name.replace('.','/'), table_file.read())
    mail.send(msg)
    for name in output_names:
        os.system('rm '+output_path+name)
    return

def features():
    session['feature1'] = ''
    session['feature2'] = ''

if __name__ == "__main__":
    app.run(port=8001)
