from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import StringField, PasswordField, BooleanField, IntegerField, TextAreaField, SubmitField, MultipleFileField,SelectMultipleField,SelectField,RadioField
from wtforms.validators import DataRequired, Length, ValidationError, Email

class geneUpload(FlaskForm):
    file = FileField('Upload File', validators=[FileRequired(), FileAllowed(['zip'])])
    predclass = StringField('Group_information', validators=[DataRequired()])
    con = StringField('Control', validators=[DataRequired()])
    test = StringField('Test Dataset File Name', validators=[DataRequired()])
    clinical = StringField('Clinical Data File Name', validators=[DataRequired()])
    quantitative = StringField('Quantitative Data Colnames', validators=[DataRequired()])
    split = StringField("interval：separated with '-'", validators=[DataRequired()])
    submit = SubmitField()

class matrixUpload(FlaskForm):
    file = FileField('Upload File', validators=[FileRequired(), FileAllowed(['csv','xlsx','txt'])])
    predclass = StringField('Group_information', validators=[DataRequired()])
    con = StringField('Control', validators=[DataRequired()])
    submit = SubmitField()



class download(FlaskForm):
    submit = SubmitField('Download')
    
class patentAssistant(FlaskForm):
    # analysis = SelectMultipleField('Based on Which Analysis',choices = [('matrix','Standard Matrix Analysis'),
    #                                                                    ('geo','Geo Dataset Analysis')])
    analysis = SelectField('Based on Which Analysis',choices = [('gene','Gene Expression Matrix'),
                                                                       ('matrix','Standard Matrix Analysis')])
    # analysis = RadioField('Based on Which Analysis',choices = [('matrix','Standard Matrix Analysis'),
    #                                                                    ('geo','Geo Dataset Analysis')],
    #                      validators=[DataRequired()])
    project = StringField('Project Name', validators=[DataRequired()])
    geneID = StringField('Gene ID', validators=[DataRequired()])
    predclass = StringField('Group_information', validators=[DataRequired()])
    con = StringField('Control', validators=[DataRequired()])
    combinationbegin = StringField('combination begin', validators=[DataRequired()])
    combinationend = StringField('combination end', validators=[DataRequired()])
    submit = SubmitField()