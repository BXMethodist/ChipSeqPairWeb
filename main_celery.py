from flask import Flask, render_template, request
from queryUtils import GEO_query
from search import Search
from match import Match
from setup import get_settings
from flask_celery import make_celery

app = Flask(__name__)
app.config['CELERY_BROKER_URL'] = 'amqp://localhost//'

celery = make_celery(app)

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

@app.route('/search',methods=['POST'])
def search():
    _searchterms = request.form['searchterms']
    _Species = request.form['Species']
    _inputEmail = request.form['inputEmail']

    settings = get_settings()
    GSMGSE_pkl = settings['GSMGSE_pkl_path']

    keywords = _searchterms.split(",")
    output_prefix = keywords[0]
    output_path = './tmp/'

    cwd = settings['Chipseq']

    CallSearch.delay(output_prefix, output_path,
               keywords, _Species, GSMGSE_pkl, cwd, _inputEmail)

    return 'We are processing your request, results will be sent to your email'

@app.route('/match',methods=['POST'])
def match():
    _feature1 = request.form['feature1']
    _feature2 = request.form['feature2']
    _Species = request.form['Species']

    _inputEmail = request.form['inputEmail']

    settings = get_settings()
    GSMGSE_pkl = settings['GSMGSE_pkl_path']

    keywords1 = _feature1.split(",")
    keywords2 = _feature2.split(",")
    output_prefix1 = keywords1[0]
    output_prefix2 = keywords2[0]
    output_path = './tmp/'

    species = _Species if _Species != '' else 'Homo sapiens'
    cwd = settings['Chipseq']

    CallMatch.delay(output_prefix1, output_prefix2, output_path,
                          keywords1, keywords2,
                          species, GSMGSE_pkl, cwd, _inputEmail)
    return 'We are processing your request, results will be sent to your email'

@app.route('/query',methods=['GET', 'POST'])
def query():
    f = request.files['IDlist']
    info = [str(line) for line in f.readlines()]
    _inputEmail = request.form['inputEmail']
    id_list = []
    for line in info:
        id_list.append(line.strip())
    f.close()
    output_path = './tmp/query.txt'
    settings = get_settings()
    GSMGSE_pkl = settings['GSMGSE_pkl_path']
    GSM_SRR_pkl = settings['GSMtoSRRpkl']

    CallQuery.delay(id_list, output_path, GSMGSE_pkl, GSM_SRR_pkl, _inputEmail)
    return 'We are processing your request, results will be sent to your email'

@celery.task(name='main_celery.search')
def CallSearch(output_prefix, output_path, keywords, _Species, GSMGSE_pkl, cwd, _inputEmail):

    Search(output_prefix, output_path, keywords, _Species, GSMGSE_pkl, cwd, _inputEmail)
    return

@celery.task(name='main_celery.match')
def CallMatch(output_prefix1, output_prefix2, output_path, keywords1, keywords2, species, GSMGSE_pkl, cwd, _inputEmail):
    Match(output_prefix1, output_prefix2, output_path, keywords1, keywords2, species, GSMGSE_pkl, cwd, _inputEmail)
    return

@celery.task(name='main_celery.query')
def CallQuery(id_list, output_path, GSMGSE_pkl, GSM_SRR_pkl, email):
    GEO_query(id_list, output_path, GSMGSE_pkl, GSM_SRR_pkl, email)
    return

if __name__ == "__main__":
    app.run()

