from flask import Flask, render_template, request, session, redirect
from queryUtils import GEO_query
from search import Search
from match import Match
from setup import get_settings
from flask_celery import make_celery
from flask_mail import Mail, Message
import json, bisect
import os

app = Flask(__name__)
app.config['CELERY_BROKER_URL'] = 'amqp://localhost//'

# app.config.update(
#     DEBUG = True,
#     MAIL_SERVER = 'ismtp.tmh.tmhs',
#     MAIL_PORT = 80,
#     MAIL_USE_SSL=True,
#     MAIL_USE_TLS=False,
#     MAIL_DEFAULT_SENDER=('Chen Lab', 'bxia@houstonmethodist.org'),
#     # MAIL_MAX_EMAIL=10,
#     MAIL_USERNAME='bxia@houstonmethodist.org',
#     MAIL_PASSWORD='20120330Xb.'
#
# )

app.config.update(
    DEBUG = False,
    MAIL_SERVER = 'smtp.gmail.com',
    MAIL_PORT = 465,
    MAIL_USE_SSL=True,
    MAIL_DEFAULT_SENDER=('Chen Lab', 'chenlabhoustonmethodist@gmail.com'),
    MAIL_MAX_EMAIL=10,
    MAIL_USERNAME='chenlabhoustonmethodist@gmail.com',
    MAIL_PASSWORD='R10-414D'
)

app.secret_key = 'R10-414D'

web_regions = open('100_10_2200_web_map.txt', 'r')
refmap = json.load(web_regions)
web_regions.close()

regions = set()
for chromosome in refmap.keys():
    # print chromosome
    for value in refmap[chromosome].values():
        regions.add((chromosome, int(value[0]), int(value[1]), int(value[2])))

refmap = {}
for region in regions:
    if region[0] in refmap:
        refmap[region[0]][region[1]] = tuple(region[1:])
        refmap[region[0]][region[2]] = tuple(region[1:])
    else:
        refmap[region[0]] = {}
        refmap[region[0]][region[1]] = tuple(region[1:])
        refmap[region[0]][region[2]] = tuple(region[1:])

link1p1 = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position='
link1p2 = '&hgsid=592956077_CVVkc6h19x28wA4sfPbMFib1ahy8'

link2p1 = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position='
link2p2 = '&hgsid=592956655_AjzBZdhCumAb8s6VcTs8mnMTyh7n'



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
    session['chr'] = 'chr1'
    session['start'] = '100'
    session['end'] = '10000'
    chr_sizes = open('chr_sizes.txt', 'r')
    session['sizes'] = json.load(chr_sizes)
    chr_sizes.close()

    peaks1 = find_regions(refmap, session['chr'], int(session['start']), int(session['end']))
    # print peaks1
    session['peaks1'] = peaks1

    session['peaks2'] = peaks1

    # print session['sizes']
    return render_template('display.html')

@app.route('/double')
def showDouble():
    print 'double'
    return render_template('double.html')

@app.route('/single')
def showSingle():
    return render_template('single.html')

@app.route('/double', methods=['POST'])
def reloadDouble():
    # print refmap.keys()
    # for key in request.form.keys():
    # print request.form
    # print request.form['L1.5x'] == 'L1.5x'
    # print request.form
    if 'submit' in request.form and request.form['submit'] == 'Submit':
        # print 'submit'
        chromosome = request.form['chr'] if request.form['chr'] != '' else session['chr']
        start = request.form['start'] if request.form['start'] != '' else session['start']
        end = request.form['end'] if request.form['end'] != '' else session['end']
        if int(start) < 1:
            start = '1'
        if int(end) > int(session['sizes'][chromosome]):
            end = str(session['sizes'][chromosome])
        session['chr'] = chromosome
        session['start'] = start
        session['end'] = end
        session['src1'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position='+chromosome+'%3A'+start+'-' + end
        session['src2'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position='+chromosome+'%3A'+start+'-' + end

    elif 'L1.5x' in request.form and request.form['L1.5x'] == 'L1.5x':
        chromosome = session['chr']
        start = session['start']
        end = session['end']

        mid = (int(end) + int(start)) / 2

        next_start = str(mid - int((int(end) - int(start))*2/3/2)/50*50)
        next_end = str(mid + int((int(end) - int(start))*2/3/2)/50*50)

        if int(next_end) - int(next_start) <= 0:
            pass
        else:
            start = next_start
            end = next_end

        if int(start) < 1:
            start = '1'
        if int(end) > int(session['sizes'][chromosome]):
            end = str(session['sizes'][chromosome])

        session['chr'] = chromosome
        session['start'] = start
        session['end'] = end
        session[
            'src1'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
        session[
            'src2'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
    elif 'L3x' in request.form and request.form['L3x'] == 'L3x':
        chromosome = session['chr']
        start = session['start']
        end = session['end']

        mid = (int(end)+int(start))/2

        next_start = str(mid - int((int(end) - int(start))/3/2)/50*50) if int(end) - int(start) != 0 else start
        next_end = str(mid + int((int(end) - int(start))/3/2)/50*50) if int(end) - int(start) != 0 else end

        if int(next_end) - int(next_start) <= 0:
            pass
        else:
            start = next_start
            end = next_end

        if int(start) < 1:
            start = '1'
        if int(end) > int(session['sizes'][chromosome]):
            end = str(session['sizes'][chromosome])

        session['chr'] = chromosome
        session['start'] = start
        session['end'] = end
        session[
            'src1'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
        session[
            'src2'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
    elif 'L10x' in request.form and request.form['L10x'] == 'L10x':
        chromosome = session['chr']
        start = session['start']
        end = session['end']

        mid = (int(end) + int(start)) / 2

        next_start = str(mid - int((int(end) - int(start))/10/2)/50*50) if int(end) - int(start) != 0 else start
        next_end = str(mid + int((int(end) - int(start))/10/2)/50*50) if int(end) - int(start) != 0 else end

        if int(next_end) - int(next_start) <= 0:
            pass
        else:
            start = next_start
            end = next_end

        if int(start) < 1:
            start = '1'
        if int(end) > int(session['sizes'][chromosome]):
            end = str(session['sizes'][chromosome])

        session['chr'] = chromosome
        session['start'] = start
        session['end'] = end
        session[
            'src1'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
        session[
            'src2'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
    elif 'S1.5x' in request.form and request.form['S1.5x'] == 'S1.5x':
        chromosome = session['chr']
        start = session['start']
        end = session['end']

        mid = (int(end) + int(start)) / 2

        next_start = str(mid - int((int(end) - int(start)) * 3 / 2 / 2) / 50 * 50)
        next_end = str(mid + int((int(end) - int(start)) * 3 / 2 / 2) / 50 * 50)

        if int(next_end) - int(next_start) <= 0:
            pass
        else:
            start = next_start
            end = next_end

        if int(start) < 1:
            start = '1'
        if int(end) > int(session['sizes'][chromosome]):
            end = str(session['sizes'][chromosome])

        session['chr'] = chromosome
        session['start'] = start
        session['end'] = end
        session[
            'src1'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
        session[
            'src2'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
    elif 'S3x' in request.form and request.form['S3x'] == 'S3x':
        chromosome = session['chr']
        start = session['start']
        end = session['end']

        mid = (int(end) + int(start)) / 2

        next_start = str(mid - int((int(end) - int(start)) * 3 / 2) / 50 * 50) if int(end) - int(start) != 0 else start
        next_end = str(mid + int((int(end) - int(start)) * 3 / 2) / 50 * 50) if int(end) - int(start) != 0 else end

        if int(next_end) - int(next_start) <= 0:
            pass
        else:
            start = next_start
            end = next_end

        if int(start) < 1:
            start = '1'
        if int(end) > int(session['sizes'][chromosome]):
            end = str(session['sizes'][chromosome])

        session['chr'] = chromosome
        session['start'] = start
        session['end'] = end
        session[
            'src1'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
        session[
            'src2'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
    elif 'S10x' in request.form and request.form['S10x'] == 'S10x':
        chromosome = session['chr']
        start = session['start']
        end = session['end']

        mid = (int(end) + int(start)) / 2

        next_start = str(mid - int((int(end) - int(start)) * 10 / 2) / 50 * 50) if int(end) - int(start) != 0 else start
        next_end = str(mid + int((int(end) - int(start)) * 10 / 2) / 50 * 50) if int(end) - int(start) != 0 else end

        if int(next_end) - int(next_start) <= 0:
            pass
        else:
            start = next_start
            end = next_end

        if int(start) < 1:
            start = '1'
        # print int(end), int(session['sizes'][chromosome])
        if int(end) > int(session['sizes'][chromosome]):
            end = str(session['sizes'][chromosome])

        session['chr'] = chromosome
        session['start'] = start
        session['end'] = end
        session[
            'src1'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
        session[
            'src2'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end

    ## current just use one map
    peaks1 = find_regions(refmap, session['chr'], int(session['start']), int(session['end']))
    # print peaks1
    session['peaks1'] = peaks1

    session['peaks2'] = peaks1

    return redirect("double")

@app.route('/single', methods=['POST'])
def reloadSingle():
    chromosome = request.form['chr']
    start = request.form['start']
    end = request.form['end']
    session['chr'] = chromosome
    session['start'] = start
    session['end'] = end
    session['src1'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position='+chromosome+'%3A'+start+'-' + end
    return render_template('reloadSingle.html')

@app.route('/display', methods=['POST'])
def selectFeatures():
    # print 'display'
    if 'feature1' in request.form or 'feature2' in request.form:
        _feature1 = request.form['feature1']
        _feature2 = request.form['feature2']

        chromosome = session['chr']
        start = session['start']
        end = session['end']

        session[
            'src1'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end
        session[
            'src2'] = 'http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=' + chromosome + '%3A' + start + '-' + end

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

@app.route('/query',methods=['GET', 'POST'])
def query():
    f = request.files['IDlist']
    info = [str(line, 'utf-8') for line in f.readlines()]
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
        table.to_csv(output_path+output_name, encoding='utf-8')
        with app.open_resource(output_path+output_name) as table_file:
            msg.attach(output_name, output_name.replace('.','/'), table_file.read())
    mail.send(msg)
    for name in output_names:
        os.system('rm '+output_path+name)
    return

def features():
    session['feature1'] = ''
    session['feature2'] = ''

def find_regions(modificationmap, chr, start, end):
    # print type(refmap['chr1'].keys()[0]),refmap['chr1'].keys()[0]
    index = sorted(modificationmap[chr].keys())

    total_width = end - start
    start_index = bisect.bisect(index, start)
    end_index = bisect.bisect(index, end)

    # print start_index, end_index
    if start_index == end_index:
        start_index -= 1

    peaks = set()
    for i in index[start_index:end_index]:
        peak = modificationmap[chr][i]
        # print peak
        peak_start = round((peak[0]-start)*1.0/total_width,2) if round((peak[0]-start)*1.0/total_width,2) >=0 else 0
        peak_end = round((peak[1]-start)*1.0/total_width,2) if round((peak[1]-start)*1.0/total_width,2) <=1 else 1
        # peak_color = 'blue' if peak[2] == 0 else 'red'
        cur_peak = (peak[0], peak[1], peak_start, peak_end-peak_start, peak[2])
        # print cur_peak
        peaks.add(cur_peak)
    if len(peaks) > 10:
        return []
    # print json.dumps(list(peaks))
    return json.dumps(list(peaks))

if __name__ == "__main__":
    app.run()

